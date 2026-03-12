#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import urllib.request as urlreq
from http.cookiejar import CookieJar
import base64    

def get_login(userz, passz, cookiez=None):
    # This is only relevant for sites that have a little browser popup to login to
    # It'll look like EXACTLY like the barebones box that pops up that is made by your browser that asks for a login
    
    if( cookiez is None ):
        cookiez = CookieJar(); # Make a cookie jar to put cookies in
    # END IF
    cookie_handler = urlreq.HTTPCookieProcessor(cookiez); # Prep a cookie handler
    
    # # Make a massword panager -> this did NOT work, left in case it's needed. web2oath would be optinal input as `web2oath=None`
    # web_pwdr = urlreq.HTTPPasswordMgrWithDefaultRealm();
    # # Add the user/pass for the site to log into to it
    # if( web2oath is None ):
    #     web_pwdr.add_password(None, web2login, userz, passz); # web2login as a req input if this is activated
    # else:
    #     web_pwdr.add_password(None, web2oath, userz, passz);
    #     # web_pwdr.add_password(None, 'https://airsl2.gesdisc.eosdis.nasa.gov/opendap', userz, passz);
    # # END IF

    # # Make an authentication handler
    # debug_handler = urlreq.HTTPHandler(debuglevel=10); # This does debug output on HTTP
    # debugS_handler = urlreq.HTTPSHandler(debuglevel=10); # This does debug output on HTTPS
    # auth_handler = urlreq.HTTPBasicAuthHandler(web_pwdr);
    # Build the URL opener
    url_opener = urlreq.build_opener(cookie_handler); # You'd add debug before this if you wanted it
    # END IF
    allurbase = base64.b64encode((userz+':'+passz).encode("utf-8")).decode("utf-8"); # Encode username and password in base64
    # Add header to the url opener
    url_opener.addheaders = [('Authorization', 'Basic '+allurbase),]; # More can be added as (tuple, sets) in this list
    
    # Set login stuff
    urlreq.install_opener(url_opener);
    
    return cookiez # If you return the cookiez and reuse them, then it'll auth faster
# END IF


