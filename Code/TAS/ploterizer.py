#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates a dictionary of plot settings to use
"""

def ploterizer(diction=None, **kwargs):
    from numpy.random import randint as np_random_randint
    import matplotlib.font_manager as fm
    import matplotlib.pyplot as plt
    
    if( diction is None ):
        diction = {'font':{}}; # Fire up a dict with a subdict 'font'
    else:
        diction['font'] = {}; # Fire up a subduct 'font' b/c the dict already exists
    # END IF
    #==============PLOT STANDARDIZATION STUFF==============
    diction['font']['size'] = {
        'small':19, #small font (tick mark font size) Default: 19, Big: 23
        'axisTick':23, #small font (tick mark font size) Default: 19, Big: 23
        'axisLabel':28, #medium font (labels) Default: 23, Big: 28
        'title':28, #big font (title only really) Default: 23, Big: 28
        'grandiose':32, #big font (title only really) Default: 26, Big: 32
        }; # standardized font sizes
    diction['font']['font'] = 'arial';
    diction['font']['weight'] = 'bold';
    
    diction['lineWidth'] = {
        'thick':4.00,
        'double plus':3.00,
        'plus':2.25,
        'regular plus':2.00,
        'regular':1.75,
        'small':1.25,
        'smaller':1,
        'smallest':0.8,
        }; #standardized line width Thick / DoublePlus / Plus / RegularPlus / Regular / Smol
    diction['linewidth'] = diction['lineWidth']; # Alias
    
    # PLOT_color = ['xkcd:purple','xkcd:green','xkcd:blue','xkcd:pink',
    #      'xkcd:brown','xkcd:brick red','xkcd:grey','xkcd:orange',
    #      'xkcd:teal','xkcd:midnight purple','xkcd:magenta','xkcd:goldenrod',
    #      'xkcd:forest green','xkcd:tan','xkcd:periwinkle','xkcd:salmon']; #color set for plotting
    diction['color'] = ['#1F77B4','#FF7F0E','#7256C1','#2CA02C','#D81B60','#656565',
                  '#0099A9','#3E0402','#480094','#FF4545','#004D40','#9A1A00',
                  '#224EB3','#249A7D','#FF45CA','#38FF32','#AF5D03','#00004E','#5F0038']; #color set for plotting [did some work at https://davidmathlogic.com/colorblind/#%231F77B4-%23FF7F0E-%237256C1-%232CA02C-%23D81B60-%23656565-%230099A9-%233E0402-%23480094-%23FF4545-%23004D40-%239A1A00-%23224EB3-%23249A7D-%23FF45CA-%2338FF32-%23AF5D03-%2300004E-%235F0038 to make it better!]
    diction['color'].extend([[valz] for valz in (np_random_randint(0,high=256,size=(600,3))/255).tolist()]); #make sure there's enough colors forever
    diction['lineStyle'] = ['-','--','-.',':']*25;
    
    
    diction['font']['fm'] = {
        'small':fm.FontProperties(family=diction['font']['font'], weight=diction['font']['weight'], size=diction['font']['size']['small']), #these are font properties, some plot stuff has this and it tells it all in one go #inconsistent
        'axisTick':fm.FontProperties(family=diction['font']['font'], weight=diction['font']['weight'], size=diction['font']['size']['axisTick']), #these are font properties, some plot stuff has this and it tells it all in one go #inconsistent
        'axisLabel':fm.FontProperties(family=diction['font']['font'], weight=diction['font']['weight'], size=diction['font']['size']['axisLabel']), #these are font properties, some plot stuff has this and it tells it all in one go #inconsistent
        'title':fm.FontProperties(family=diction['font']['font'], weight=diction['font']['weight'], size=diction['font']['size']['title']), #these are font properties, some plot stuff has this and it tells it all in one go #inconsistent
        'grandiose':fm.FontProperties(family=diction['font']['font'], weight=diction['font']['weight'], size=diction['font']['size']['grandiose']), #these are font properties, some plot stuff has this and it tells it all in one go #inconsistent
        }; # standardize font managers
        
    diction['width_1C'] = 3.5; #plot width for a single column plot
    diction['width_2C'] = 7.5; #plot width for a double column plot (full page width)
    diction['height_max'] = 9.8; #plot height limit (can be less)
    diction['dpi'] = 300; #pixels per inch
    diction['fileType'] = '.png'; #file type to save figure as, Default: .png, Alts: .pdf, .jpg, more...
    diction['unitBracket_L'] = '['; #whatevs you want for left bracket around units
    diction['unitBracket_R'] = ']'; #whatevs you want for right bracket around units
    
    # Alias
    diction['font']['size']['smol'] = diction['font']['size']['small']; # For backwards compatibility
    diction['font']['fm']['smol'] = diction['font']['fm']['small']; # For backwards compatibility
    diction['lineWidth']['smol'] = diction['lineWidth']['small']; # For backwards compatibility
    diction['lineWidth']['smoller'] = diction['lineWidth']['smaller']; # For backwards compatibility
    diction['lineWidth']['smollerest'] = diction['lineWidth']['smallest']; # For backwards compatibility
    diction['ext'] = diction['fileType']; # For backwards compatibility
    diction['save file type'] = diction['fileType']; # For backwards compatibility+
    
    # Override kwargs as requested
    def kwargerizer(diction, kwargs):
        for kwargz in kwargs:
            # If statement for dealing with dicts-in-dicts
            if( isinstance(kwargs[kwargz], dict) and (kwargz in diction) ):
                diction = kwargerizer(diction, kwargs); # Go deeper
            else:
                diction[kwargz] = kwargs[kwargz]; # Set it
            # END IF
        # END FOR kwargz
        return diction
    # END DEF
    
    diction = kwargerizer(diction, kwargs);
    
    # Apply defaults to help stuff be standardized easier
    plt.rcParams['font.weight'] = diction['font']['weight']; #sents default font weight to bold for everything else
    plt.rcParams['axes.labelweight'] = diction['font']['weight']; #sets default font weight to bold for axis labels
    plt.rc('font', size=diction['font']['size']['axisTick']); #default text size
    #plt.rc('axes', titlesize=diction['font']['size']['axisTick']); # axes title text size
    #plt.rc('axes', labelsize=diction['font']['size']['axisTick']); #x and y label font size
    plt.rc('xtick', labelsize=diction['font']['size']['axisTick']); #x tick label font size
    plt.rc('ytick', labelsize=diction['font']['size']['axisTick']); #y tick label font size
    plt.rc('legend', fontsize=diction['font']['size']['axisLabel']); #legend fornt size
    plt.rc('figure', titlesize=diction['font']['size']['title']); #figure title font size (this one didn't do anything, so 2nd here also)
    plt.rc('axes', titlesize=diction['font']['size']['title']); #figure title font size (this one did work)
    plt.rcParams['axes.labelsize'] = diction['font']['size']['axisLabel']; #try this one also
    #plt.rc('text', fontsize=diction['font']['size']['axisTick']); #figure text font size
    plt.rcParams['axes.unicode_minus'] = False; #disable big - sign and instead use smol - sign (saves plot space)

    return diction
# END DEF

# Class version, complexity not needed
# class ploterizer:
#     def __init__(self, **kwargs):
#         from numpy.random import randint as np_random_randint
#         import matplotlib.font_manager as fm
        
#         print(kwargs)
        
#         #==============PLOT STANDARDIZATION STUFF==============
#         self.FONT_smol = 19; #small font (tick mark font size) Default: 19, Big: 23
#         self.FONT_axisTick = 23; #small font (tick mark font size) Default: 19, Big: 23
#         self.FONT_axisLabel = 28; #medium font (labels) Default: 23, Big: 28
#         self.FONT_title = 28; #big font (title only really) Default: 23, Big: 28
#         self.FONT_grandiose = 32; #big font (title only really) Default: 26, Big: 32
#         self.FONT_font = 'arial';
#         self.FONT_weight = 'bold';
        
#         self.PLOT_lineWidth = {
#             'thick':4.00,
#             'double plus':3.00,
#             'plus':2.25,
#             'regular plus':2.00,
#             'regular':1.75,
#             'smol':1.25,
#             'smoller':1,
#             'smollerest':0.8,
#             }; #standardized line width Thick / DoublePlus / Plus / RegularPlus / Regular / Smol
#         # PLOT_color = ['xkcd:purple','xkcd:green','xkcd:blue','xkcd:pink',
#         #      'xkcd:brown','xkcd:brick red','xkcd:grey','xkcd:orange',
#         #      'xkcd:teal','xkcd:midnight purple','xkcd:magenta','xkcd:goldenrod',
#         #      'xkcd:forest green','xkcd:tan','xkcd:periwinkle','xkcd:salmon']; #color set for plotting
#         self.PLOT_color = ['#1F77B4','#FF7F0E','#7256C1','#2CA02C','#D81B60','#656565',
#                       '#0099A9','#3E0402','#480094','#FF4545','#004D40','#9A1A00',
#                       '#224EB3','#249A7D','#FF45CA','#38FF32','#AF5D03','#00004E','#5F0038']; #color set for plotting [did some work at https://davidmathlogic.com/colorblind/#%231F77B4-%23FF7F0E-%237256C1-%232CA02C-%23D81B60-%23656565-%230099A9-%233E0402-%23480094-%23FF4545-%23004D40-%239A1A00-%23224EB3-%23249A7D-%23FF45CA-%2338FF32-%23AF5D03-%2300004E-%235F0038 to make it better!]
#         self.PLOT_color.extend((np_random_randint(0,high=256,size=(200,3))/255).tolist()); #make sure there's enough colors forever
#         self.PLOT_lineStyle = ['-','--','-.',':']*25;
        
        
#         self.FONT_smolFM = fm.FontProperties(family=self.FONT_font, weight=self.FONT_weight, size=self.FONT_smol); #these are font properties, some plot stuff has this and it tells it all in one go #inconsistent
#         self.FONT_axisTickFM = fm.FontProperties(family=self.FONT_font, weight=self.FONT_weight, size=self.FONT_axisTick); #these are font properties, some plot stuff has this and it tells it all in one go #inconsistent
#         self.FONT_axisLabelFM = fm.FontProperties(family=self.FONT_font, weight=self.FONT_weight, size=self.FONT_axisLabel); #these are font properties, some plot stuff has this and it tells it all in one go #inconsistent
#         self.FONT_titleFM = fm.FontProperties(family=self.FONT_font, weight=self.FONT_weight, size=self.FONT_title); #these are font properties, some plot stuff has this and it tells it all in one go #inconsistent
#         self.FONT_grandioseFM = fm.FontProperties(family=self.FONT_font, weight=self.FONT_weight, size=self.FONT_grandiose); #these are font properties, some plot stuff has this and it tells it all in one go #inconsistent
    
#         self.journal_width_1C = 3.5; #plot width for a single column plot
#         self.journal_width_2C = 7.5; #plot width for a double column plot (full page width)
#         self.journal_height_max = 9.8; #plot height limit (can be less)
#         self.journal_dpi = 300; #pixels per inch
#         self.plot_fileType = '.png'; #file type to save figure as, Default: .png, Alts: .pdf, .jpg, more...
#         self.plot_unitBracket_L = '['; #whatevs you want for left bracket around units
#         self.plot_unitBracket_R = ']'; #whatevs you want for right bracket around units
        
#         # Override kwargs as requested
#         for kwargz in kwargs:
#             self.kwargz = kwargs[kwargz]
#         # END FOR kwargz
#     # END DEF
# # END CLASS