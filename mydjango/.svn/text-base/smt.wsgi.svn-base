import os
import sys

smthome = '/priv/maipenrai/skymap/skymap/subpipe'

if smthome not in sys.path:
    sys.path.append(smthome)

os.environ.setdefault("SUBPIPEHOME", smthome)

if os.path.isdir('/ramdisk'):
    os.environ.setdefault("SUBSCRATCH",'/ramdisk/subpipe_production')
else:
    os.environ.setdefault("SUBSCRATCH",'%s/subpipe_production'%smthome)

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "mydjango.settings")

import django.core.handlers.wsgi
application = django.core.handlers.wsgi.WSGIHandler()


