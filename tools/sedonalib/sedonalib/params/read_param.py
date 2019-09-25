from SedonaParam import *

def new_paramfile(template=None):

    param = SedonaParam()
    if (template is not None):
        param.set_template(template)
    return param
