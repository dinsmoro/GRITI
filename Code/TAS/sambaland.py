#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sambaland, like Shondaland

Version
1.0 - initial feature set
1.1 - added initial support for a path to a file

==Initialize like this in Python:
 > from sambaland import sambaland
 > smb = sambaland('//samba/share/to/mount', FLG_crashOnCollide=False)
++ sambaland takes 'smb://samba/share/to/mount', '//samba/share/to/mount', 
++     '/samba/share/to/mount' (don't do this on the reg, that's local), and 
++     Windows-style '\\samba\share\to\mount'
==Access the files in the mount like:
 > with open( os.path.join(smb.local, 'textfiletoread.txt'), r):
++You can see the files yourself by getting the smb.local string, 
++    copying it, going to Finder (in Mac), Go>Go to Folder... and 
++    pasting it the string in and hitting enter.
++The entire samba server (smb.server name without the smb://) 
++    is also under the "Locations" zone in the Mac Finder's 
++    left hand side column of links
==Close the connection with:
 > smb.close()
 
 
==Required:
    sambaShare - the path to a network Samba share in the form of
        'smb://samba/share/to/mount', '//samba/share/to/mount', 
        '/samba/share/to/mount' (don't do this on the reg), and 
        Windows-style '\\samba\share\to\mount'
==Options:
    FLG_crashOnCollide - True or False, True crashes if it looks like the the 
        samba share is already mounted, False lets it slide 
        (usually works out that folder names don't collide - but they could!)                                                                                                      
    silent - True or False, True surpresses anything below ERROR warnings 
        (which will crash the program anyway), False lets warnings chat
"""

class sambaland:
    def __init__(self, sambaShare, FLG_crashOnCollide=False, silent=False):
        import re
        from os import path as ospath#, system as ossystem
        from subprocess import run as subrun, DEVNULL
        from platform import system as platform_system
        # Get Python class busywork sorted
        self.share = sambaShare;
        self.FLG_crashOnCollide = FLG_crashOnCollide;
        self.silent = silent;
        self.live = False; # Default to not live connection
        
        # Get the OS
        self.osname = platform_system();
        
        # Initialize samba connection
        if( self.osname == 'Windows' ):
            # from warnings import warn
            # warn('Warning in sambaland.__init__: Windows is unsupported for automatically loading Samba shares, sorry.');
            if( (self.share[0:2] == '//') or (self.share[0:4] == r'smb:') ):
                self.share = self.share.replace('/', '\\').replace('smb:',''); # Convert Unix path to Windows
                if( self.share[0:2] != r'\\' ):
                    self.share = '\\\\'+self.share[0:2].replace(r'\\','')+self.share[2:];
                # END IF
            # END IF
            
            # Standardize the server connection using one line if statements
            if( self.share[-1] == '\\' ): self.share = self.share[:-1]; # Ending \\ might cause problems with the local path flipping stuff, so ditch it            
            # Windows doesn't seem to need to explicitly mount Samba shares if you've already done it; you can just pathe to them
            self.server = self.share[:re.match(r'^\\{1,2}[^\\]+\\', self.share).end()-1]; # Server location
            self.path = self.share[re.match(r'^\\{1,2}[^\\]+\\', self.share).end():]; # Path on server for the data
            self.local = self.share; # Pass through, can use directly
            self.live = True; # Note that the connection is live (Might not need since Windows doesn't have explicit mount/umount view of Samba shares)
            
        elif( self.osname == 'Darwin' ):
            # Convert to Unix-style path if it's a Windows-style path
            
            if( self.share[0:2] == r'\\' ):
                self.share = self.share.replace('\\', '/'); # Convert Windows path to Unix
            # END IF
            
            # Standardize the server connection using one line if statements
            if( re.match(r'^\/[^\/]+?\/', self.share) is not None ): self.share = 'smb:/'+self.share;
            if( re.match(r'^\/\/[^\/]+?\/', self.share) is not None ): self.share = 'smb:'+self.share;
            if( re.match(r'^smb:[^\/]+?\/', self.share) is not None ): self.share = 'smb://'+self.share[4:];
            if( self.share[-1] == '/' ): self.share = self.share[:-1]; # Ending / causes problems with the local path flipping stuff, so ditch it
            
            if( ospath.splitext(self.share)[1] == '' ): # If a possible extension is empty, go ahead and guess
                # Discern relevant locations/sublocations
                self.server = self.share[:re.match(r'(?:^|^smb:)\/{1,2}[^\/]+?\/', self.share).end()-1]; # Server location
                self.path = self.share[re.match(r'(?:^|^smb:)\/{1,2}[^\/]+?\/', self.share).end():]; # Path on server for the data
                self.local = ospath.join('/Volumes', self.path[::-1][:re.match(r'^[^\/]+?\/', self.path[::-1]).end()-1][::-1]); # Where the data will be locally mounted on MacOS
                
                if( not ospath.isdir(self.local) ):
                    # Now mount the samba share
                    #ossystem('osascript -e \'mount volume "'+self.share+'"\' >/dev/null 2>1');
                    subrun('osascript -e \'mount volume "'+self.share+'"\'', shell=True, stdout=DEVNULL, check=True); # Mount the samba share to the local computer using the built-in Mac stuff (gets to use the saved user/password if you've connected to it once before!)
                    self.live = True; # Note that the connection is live
                else:
                    if( self.FLG_crashOnCollide == True ):
                        raise Exception('ERROR in sammbaland: Samba share `'+self.share+'` seems to be already mounted, or worse there is a name collision. `FLG_crashOnCollide` is set to `True` so raising an error.');
                    else:
                        self.live = True; # Assume the connection is OK and live
                        if( self.silent == False ):
                            print('WARNING in sammbaland: Samba share `'+self.share+'` seems to be already mounted, or worse there is a name collision. `FLG_crashOnCollide` is not set to `True` so continuing on.');
                        # END IF
                    # END IF
                # END IF
            else:
                # If there is a possible extension, it may be a file. Try to fall back one level to guarantee not-a-file
                self.server = self.share[:re.match(r'(?:^|^smb:)\/{1,2}[^\/]+?\/', self.share).end()-1]; # Server location
                self.path = self.share[re.match(r'(?:^|^smb:)\/{1,2}[^\/]+?\/', self.share).end():]; # Path on server for the data
                self.local = ospath.join('/Volumes', ospath.dirname(self.share)[::-1][:re.match(r'^[^\/]+?\/', ospath.dirname(self.share)[::-1]).end()-1][::-1], self.path[::-1][:re.match(r'^[^\/]+?\/', self.path[::-1]).end()-1][::-1]); # Where the data will be locally mounted on MacOS
                
                if( not ospath.isdir(ospath.join('/Volumes', ospath.dirname(self.share)[::-1][:re.match(r'^[^\/]+?\/', ospath.dirname(self.share)[::-1]).end()-1][::-1])) ): # Heuristic to avoid remounting assuming it is a file path
                    subrun('osascript -e \'mount volume "'+ospath.dirname(self.share)+'"\'', shell=True, stdout=DEVNULL, check=True); # Mount the samba share to the local computer using the built-in Mac stuff (gets to use the saved user/password if you've connected to it once before!)
                    self.live = True; # Note that the connection is live
                else:
                    if( self.FLG_crashOnCollide == True ):
                        raise Exception('ERROR in sammbaland: Samba share `'+self.share+'` seems to be already mounted, or worse there is a name collision. `FLG_crashOnCollide` is set to `True` so raising an error.');
                    else:
                        self.live = True; # Assume the connection is OK and live
                        if( self.silent == False ):
                            print('WARNING in sammbaland: Samba share `'+self.share+'` seems to be already mounted, or worse there is a name collision. `FLG_crashOnCollide` is not set to `True` so continuing on.');
                        # END IF
                    # END IF
                # END IF
            # END IF
        elif( self.osname == 'Linux' ):
            # HPC doesn't use Samba directly, so no way to test this atm
            # from warnings import warn
            # warn('Warning in sambaland.__init__: Linux is unsupported for automatically loading Samba shares, sorry.');
            
            if( self.share[0:2] == r'\\' ):
                self.share = self.share.replace('\\', '/'); # Convert Windows path to Unix
            # END IF
            
            self.local = self.share; # Pass through for now
            self.path = self.share;
            self.server = self.share;
        else:
            from warnings import warn
            warn('Warning in sambaland.__init__: OS name detected is unrecognized. Printing detected OS name: '+self.osname);
            self.local = self.share; # Pass through for now
            self.path = self.share;
            self.server = self.share;
        # END IF
    # END DEF
    
    def close(self):
        if( self.live == True ):
            if( self.osname == 'Windows' ):
                from warning import warn
                warn('Warning in sambaland.close: Closing a samba connection on Windows is not supported right now, sorry.');
            elif( self.osname == 'Darwin' ):
                #from os import system as ossystem
                #ossystem('umount '+self.local+' >/dev/null 2>&1')
                from subprocess import run as subrun, DEVNULL
                subrun('umount '+self.local, shell=True, stdout=DEVNULL, check=True); # Unmount the samba share volume from the local computer
                self.live = False; # Record connection is donezo
            elif( self.osname == 'Linux' ):
                from warning import warn
                warn('Warning in sambaland.close: Closing a samba connection on Linux is not supported right now, sorry.');
            else:
                from warnings import warn
                from warning import warn
                warn('Warning in sambaland.close: OS name detected is unrecognized. Printing detected OS name: '+self.osname);
            # END IF
            # END IF
        else:
            if( self.silent == False ):
                print('WARNING in sammbaland: Samba share `'+self.share+'` has already been closed.');
            # END IF
        # END IF
    # END DEF
# END CLASS

# smbprotocol-based solution (not needed on Mac, may want as backup(?))
# note smbprotocol needs gssapi and krb5 as well
# import smbclient
# from getpass import getuser, getpass
# from os import environ as os_environ
# import warnings

# # --- Samba Support ---
# def samba_credentials(smb_pass=None):
#     if 'SPY_PYTHONPATH' in os_environ:
#         def prompt_password(user, reason='Credential'):
#             """
#             Parameters
#             ----------
#             user : user name
        
#             Returns
#             -------
#             text : user input password
#             """
#             # Guaranteed installed because Spyder needs it
#             from PyQt5.QtWidgets import  QInputDialog, QLineEdit, QApplication, QWidget
#             from PyQt5.QtCore import Qt, QCoreApplication
        
#             # Let's avoid to crash Qt-based dev environment (Spyder...)
#             app = QCoreApplication.instance();
#             if app is None:
#                 app = QApplication([]);
#             # END IF
            
#             FLG_donezo = False;
#             while( FLG_donezo == False ):
#                 text, ok = QInputDialog.getText(None, reason, 'Enter password for {}:'.format(user), QLineEdit.Password, flags=Qt.WindowStaysOnTopHint);
#                 if ok and text:
#                     FLG_donezo = True;
#                     return text
#                 else:
#                     warnings.warn('WARNING: You must enter a password and it must be text.')
#                 # END IF
#             # END IF
#             #raise ValueError('Must enter a password and it must be text.');
#         # END DEF
#     # END IF
    
#     # Get user name
#     smb_user = getuser()
#     print(smb_pass)
#     # Ask password if not already present from re-runs
#     if( smb_pass != True ):
#         if( 'SPY_PYTHONPATH' in os_environ ):
#             smb_pass = prompt_password(smb_user, 'Samba Password');
#         else:
#             # Spyder cannot do hidden text
#             print('Enter password for user '+smb_user+':')
#             smb_pass = getpass();
#         # END IF
        
#         # great reference: https://scriptingosx.com/2018/08/user-interaction-from-bash-scripts/
#         #smb_pass = subprocess.getoutput('osascript -e \'display dialog "Enter password for '+smb_user+':" with hidden answer\'');
#         #smb_pass = subprocess.run(['osascript -e \'display dialog "Enter password for '+smb_user+':" with hidden answer\'', capture_output=True, text=True);
#         #smb_pass = subprocess.check_output(['osascript', '-e \'display dialog "Enter password for '+smb_user+':" with hidden answer\'']);
#         # subprocess.getoutput("""osascript -e 'display dialog "Who are you?" default answer "nobody" with hidden answer'""");
        
#         # Apply
#         smbclient.ClientConfig(username=smb_user, password=smb_pass);
        
#         # Remove plaintext password
#         del smb_pass
        
#         return True
#     # END IF
    # # END DEF