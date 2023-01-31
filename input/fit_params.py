import gvar as gv

fit_params = {
    'S_12':{# S_12 is also known as D(p,gamma)3-He
        'dsets'     :['Casella_2002', 'Ma_1997', 'Mossa_2020', 'Schmid_1997',
                      'Tisma_2019', 'Turkat_2021', 'Warren_1963'],
        'plot_params':{
            'x_lim'     :(0.00065, 2.5),
            'y_lim'     :(8e-08, 2.5e-05),
            'x_scale'   :'log',
            'y_scale'   :'log',
            'x_label'   :'Energy (MeV)',
            'y_label'   :'S-factor (MeV b)',
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
        'models':['poly 2', 'poly 3', 'poly 4', 'poly 5'],
        'priors':{
            'S_0':gv.gvar(2e-7,1e-7),
            'S_1':gv.gvar(0, 1),
            'S_2':gv.gvar(0, 1),
            'S_3':gv.gvar(0, 1),
            'S_4':gv.gvar(0, 1),
            'S_5':gv.gvar(0, 1)},
    }
}
