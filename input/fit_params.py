import gvar as gv
import numpy as np

fit_params = {
    'S_12':{# S_12 is also known as D(p,gamma)3-He
        'dsets'      :['Casella_2002', 'Ma_1997', 'Mossa_2020', 'Schmid_1997',
                      'Tisma_2019', 'Turkat_2021', 'Warren_1963'],
        'plot_params':{
            'x_lim'     :(0.00065, 2.5),
            'y_lim'     :(8e-08, 2.5e-05),
            'x_scale'   :'log',
            'y_scale'   :'log',
            'x_label'   :'Energy (MeV)',
            'y_label'   :'S factor (MeV b)',
            'text'      :r'D(p,$\gamma$)$^3$He',
            'text_loc'  :[0.1,0.6],
            'legend_loc':4,
            'colors'    :{
                'Casella_2002':'blueviolet',
                'Ma_1997'     :'saddlebrown',
                'Mossa_2020'  :'r',
                'Schmid_1997' :'orange',
                'Tisma_2019'  :'b',
                'Turkat_2021' :'silver',
                'Warren_1963' :'g',
            },
        },
        'models':['pheno', 'pheno_offset', 'poly_2', 'poly_3', 'poly_4', 'poly_5'],
        'priors':gv.BufferDict({
            'S_0':gv.gvar(2e-7,1e-7),
            'S_1':gv.gvar(0, 1e-5),
            'S_2':gv.gvar(0, 1e-5),
            'S_3':gv.gvar(0, 4e-6),
            'S_4':gv.gvar(0, 4e-6),
            'S_5':gv.gvar(0, 4e-6),
            'S_6':gv.gvar(0, 4e-6),
            'a'  :gv.gvar(1, 0.08),
            'b'  :gv.gvar(0, 1e-8),
            'log(f_Turkat_2021)' :gv.gvar(0.00, np.log(1.14)),
            'log(f_Mossa_2020)'  :gv.gvar(0.00, np.log(1.027)),
            'log(f_Tisma_2019)'  :gv.gvar(0.00, np.log(1.10)),
            'log(f_Casella_2002)':gv.gvar(0.00, np.log(1.045)),
            'log(f_Schmid_1997)' :gv.gvar(0.00, np.log(1.09)),
            'log(f_Ma_1997)'     :gv.gvar(0.00, np.log(1.09)),
            'log(f_Warren_1963)' :gv.gvar(0.00, np.log(1.10))
            }),
        'z0'    :{'Turkat_2021':0.1, 'Mossa_2020':0.1,
                    'Tisma_2019':0.1, 'Casella_2002':0.1, 'Schmid_1997':0.1,
                    'Ma_1997':0.1,'Warren_1963':0.1}
    }
}
