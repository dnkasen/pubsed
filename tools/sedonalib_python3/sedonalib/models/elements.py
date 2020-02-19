import numpy as np

def get_composition(name,elems=None,Xlan=None):


    if (name == "solar"):

        this_elem = ["1.1" ,"2.4"]
        this_comp = [0.75,   0.25]

#    if (name == "r-process"):



def get_rprocess_composition(Xf=None,Xd = None):

    # default f-shell abundance is 10%
    if (Xf is None):
        Xf = 0.1

    if (Xd is None):
        Xd = 1 - Xf

    Xs = 1.0 - Xf - Xd
    Xp = 0.0

    Z_f = np.array([58,59,60,61,62,63,64,65,66,67,68,69,70])
    A_f = np.array([140,140,144,144,150,151,157,158,162,168,168,173,173])
    C_f = np.array([0.00275810142673,  0.00114986769631,  0.0201991683927,  0.002,  0.00645021110853,  0.00243513357583,  0.00857752936147,  0.00154915757877,  0.0103926469702,  0.00241689087521,  0.00707442508677,  0.00103760326605,  0.00711672655612])
    C_f = C_f/sum(C_f)

    Z_d = np.array([21,22,23,24, 25, 26, 27, 28],dtype='i')
    A_d = np.array([92,94,96,97,100,102,104,104],dtype='i')
    C_d = np.array([ 1.,1.,1.,1., 1.,1., 1., 1.])
    C_d = C_d/sum(C_d)

    Z_p = np.array([])
    A_p = np.array([])
    C_p = np.array([])
    C_p = C_p/sum(C_p)

    Z_s = np.array([20,38])
    A_s = np.array([90,90])
    C_s = np.array([1.0,0.00])
    C_s = C_s/sum(C_s)

    if (Xd == 0):
        Z_d = []
        A_d = []
        C_d = []
    if (Xs == 0):
        Z_s = []
        A_s = []
        C_s = []
    if (Xf == 0):
        Z_f = []
        A_f = []
        C_s = []

    Z = np.concatenate((Z_s,Z_p,Z_d,Z_f))
    A = np.concatenate((A_s,A_p,A_d,A_f))
    comp = np.concatenate((C_s*Xs,C_p*Xp,C_d*Xd,C_f*Xf))
    comp = comp/sum(comp)

    this_elem = []
    for this_Z,this_A in zip(Z,A):
        this_elem.append("{0:.0f}.{1:.0f}".format(this_Z,this_A))

    return this_elem,comp


def get_element_ZA(e):

    if (not isinstance(e, str)):
        e = str(e)


    # read Z.A format
    if ("." in e):
        Z,A =  e.split(".")

    # read name format
    else:

        Z_dict = {
            "H":  1,
            "He": 2,
            "Li": 3,
            "C": 6,
            "O": 8
        }

        temp = re.findall(r'\d+', e)
        A = list(map(int, temp))
        if (len(A) !=1):
            raise exception("Error: " + e + " is not a valid element name\n")
        A = A[0]
        ename = e.replace(str(A),"")
        Z = Z_dict[ename]

    A = int(A)
    Z = int(Z)
    return Z,A