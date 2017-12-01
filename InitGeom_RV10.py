import CalcAircraft.Geometry as geom
import CalcAircraft.CEA_VLM as CEA_VLM
import CalcAircraft.Stability as stab
import numpy as np
import pdb
import math
import scipy.interpolate as interpolate
import scipy.integrate as integrate

def load_geom(x = []):
    

    
    ### Inicializa aviao
    airc_data = dict()
    airc = geom.Geom()
    ######## Dados relativos a asa #########
    
    # Carrega pontos do BA e BF de arquivo de texto
    
    print('Carregando dados da asa...\n')
    BA_tmp = np.loadtxt('Geometry_files/RV10_LE.txt')
    BF_tmp = np.loadtxt('Geometry_files/RV10_TE.txt')
    
    BA_tmp2 = np.flipud(BA_tmp) * [1,-1,1]
    BF_tmp2 = np.flipud(BF_tmp) * [1,-1,1]
    airc.wings[0].BA = np.array(list(BA_tmp2[:len(BA_tmp2)-1,:])+list(BA_tmp))
    airc.wings[0].BF = np.array(list(BF_tmp2[:len(BA_tmp2)-1,:])+list(BF_tmp))
    airc.wings[0].calc_area_CMA()
    
    wing = dict()
    wing['cr'] = 1.420 #corda na raiz
    wing['pos'] = 2.05345 #Distancia bordo de ataque da corda da raiz ao datum
    wing['incidencia'] = 1.0 #Incidência da Asa
    # ht['wing_altura'] = -0.4143370 #Altura da Asa em Relacao ao Datum 'Spinner'
    
    wing['airf_r'] = 'RV10.dat' 
    wing['airf_afil'] = 'RV10.dat' 
    wing['airf_trans'] = 'RV10.dat'
    wing['airf_tip'] = 'RV10.dat' 
    
    #Definindo posicoes pelo length
    centro = round(len(airc.wings[0].BF)/2)
    i_fim = len(airc.wings[0].sec_L) - 1
    wing['L_r'] = airc.wings[0].sec_L[centro+1]
    wing['L_afil'] = airc.wings[0].sec_L[i_fim-2]
    wing['L_trans'] = airc.wings[0].sec_L[i_fim-1]
    wing['L_tip'] = np.max(airc.wings[0].sec_L)
    
    perfis =  [wing['airf_tip'],wing['airf_trans'],wing['airf_afil'],wing['airf_r'],wing['airf_r'],wing['airf_afil'],wing['airf_trans'],wing['airf_tip']]
    posicao = [-wing['L_tip'],-wing['L_trans'],-wing['L_afil'],-wing['L_r'],wing['L_r'],wing['L_afil'],wing['L_trans'],wing['L_tip']]
    
    # Atribuicao de perfis nas posicoes calculadas
    airc.wings[0].select_airf(posicao,perfis,calc=True, length = True,symmetry = False)
    airc.wings[0].n_pan([-np.max(airc.wings[0].sec_L), np.max(airc.wings[0].sec_L)],[20],length = True,symmetry = False)
    
    
    ## Carrega dados da empenagem horizontal
    ht = dict()
    ht['pos'] = 7.0356502    # Distancia do BA ao datum
    ht['altura'] = 0.7776459 # Altura do BA em relaçao ao BA da Asa
    
    print('Carregando dados da EH...\n')
    BAh_tmp = np.loadtxt('Geometry_files/RV10_EH_LE.txt')
    BFh_tmp = np.loadtxt('Geometry_files/RV10_EH_TE.txt')
    
    BAh_tmp2 = np.flipud(BAh_tmp) * [1,-1,1]
    BFh_tmp2 = np.flipud(BFh_tmp) * [1,-1,1]
    airc.wings[1].BA = np.array(list(BAh_tmp2[:len(BAh_tmp2)-1,:])+list(BAh_tmp))
    airc.wings[1].BF = np.array(list(BFh_tmp2[:len(BAh_tmp2)-1,:])+list(BFh_tmp))
    airc.wings[1].calc_area_CMA()
    
    
    perfis_ht = ['NACA0012.dat', 'NACA0012.dat','NACA0012.dat','NACA0012.dat'] #RV10 utiliza o NACA0012
    posicao_ht = [airc.wings[1].BA[0,1],airc.wings[1].BA[1,1],airc.wings[1].BA[2,1],airc.wings[1].BA[3,1]]
    
    airc.wings[1].select_airf(posicao_ht,perfis_ht,symmetry=False,calc=True)
    airc.wings[1].new_flap([-max(airc.wings[1].BA[:,1]) + 0.01,-0.05],[0.50,0.50],b_trans = [0,0],deflex = 0)
    airc.wings[1].new_flap([0.05,max(airc.wings[1].BA[:,1]) - 0.01],[0.50,0.50],b_trans = [0,0],deflex = 0)
    airc.wings[1].n_pan([-np.max(airc.wings[1].sec_L), np.max(airc.wings[1].sec_L)],[10],length = True,symmetry = False)
    print('Polares 2D calculadas da EH.\n')
    airc.lon_control_surf = [1]
    airc.wings[1].lon_control = list(range(len(airc.wings[1].flaps)))
    
    
    ## Carrega dados da fuselagem
    fus = dict()
    print('Carregando dados da Fuselagem...\n')
    fus['D_w'] = 0.5847758*2 #Diâmetro da Fuselagem na Secao do Bordo de Ataque da Asa
    fus['D_ht'] = 0.1663523*2 #Diâmetro da Fuselagem na Secao do Bordo de Ataque da Empenagem H. 
    fus['alfa0'] = -3 #Incidencia da Fuselagem em relacao ao Angulo de Sustentacao Nula da Fuselagem (Bom senso segundo Pullin)
    fus['rug'] = 6.4e-6  #Rugosidade Manter Padrao
    x_fus = np.array([0,0.3769625,0.5084285,0.7206020,1.2224689,2.0468630,2.6597150,3.7167103,4.8437909,6.0852766,7.0356502,7.4422]) #Coordenadas em X das Secoes ao Longo da Fuselagem em relacao ao Datum
    larg_fus = np.array([0,0.1919531,0.4387350,0.4733639,0.5255099,0.5847758,0.6116011,0.5709773,0.3910283,0.1663523,0.0362937,0]) #Largura das Secoes ao Longo da Fuselagem
    p_fus = np.array([0,1.925,2.373,3.113,3.502,4.110,4.425,4.006,3.128,1.857,1.040,0]) #Perimetro das Secoes ao Longo da Fuselagem
    fun_larg = interpolate.interp1d(x_fus,larg_fus,kind='cubic')
    fun_p = interpolate.interp1d(x_fus,p_fus,kind='cubic')
    fus['lb'] = 4.7824850 #Comprimento da Secao B em Relacao ao termino da Secao A
    fus['lA'] = 2.6539214 #Comprimento da Secao A em Relacao ao Datum 
    fus['Sb'] = integrate.quad(fun_larg,min(x_fus),fus['lb'])[0] 
    fus['enflex_bf'] = 3 #Angulo entre o eixo de referencia da fuselagem e a regiao da cauda (pag 64)
    fus['lbf'] = wing['pos']+wing['cr']/4.0
    fus['Sbf'] = integrate.quad(fun_larg,min(x_fus),fus['lbf'])[0]
    fus['Sm'] = integrate.quad(fun_p,min(x_fus),fus['lb'])[0]
    fus['Smc'] = integrate.quad(fun_p,fus['lA'],fus['lb'])[0]
    fus['Scab'] = 0.437 #Area Frontal do Parabrisa 
    fus['x0'] = (0.374+0.533*fus['lA']/fus['lb'])*fus['lb']
    fus['Sb_x0'] = (integrate.quad(fun_p,fus['x0'],fus['lb'])[0])/np.pi
        
    # Dados gerais da aeronave
    airc_data['h'] = 0.24
    airc.CG = [airc_data['h']*airc.wings[0].CMA,0,0]
    
    airc.wing_data = wing
    airc.fus_data  = fus
    airc.ht_data   = ht
    airc.airc_data = airc_data
    
    airc.wings[0].plot_3D() #######
    airc.wings[1].plot_3D() #######


    
    return airc
