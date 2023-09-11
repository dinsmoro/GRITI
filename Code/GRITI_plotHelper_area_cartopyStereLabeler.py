#labels cause cartopy doesn't do it for some reason
import cartopy.crs as ccrs #cartopy replaces basemap b/c it's updated
from Code.subfun_textNice import textNice

def GRITI_plotHelper_area_cartopyStereLabeler(ax,ticks_lat,ticks_long,settings_plot,min_lat,FLG_latLabel=True,FLG_longLabel=True):
    if( FLG_longLabel == True ):
        for j in range(0,ticks_long.size):
            ax.text(ticks_long[j], min_lat, textNice(ticks_long[j])+'°', fontproperties=settings_plot['font axis tick FM'], \
                       horizontalalignment='center', transform=ccrs.PlateCarree(), clip_on=False, zorder=375);
        #END FOR j
    #END IF
    if( FLG_latLabel == True ):
        for j in range(0,ticks_lat.size):
            ax.text(135, ticks_lat[j], textNice(ticks_lat[j])+'°', fontproperties=settings_plot['font axis tick FM'], \
                       horizontalalignment='center', transform=ccrs.PlateCarree(), clip_on=False, zorder=375);
        #END FOR j
    #END IF
#END DEF