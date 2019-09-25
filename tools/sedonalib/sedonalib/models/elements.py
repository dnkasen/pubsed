

def get_composition(name):

    if (name == "solar"):
        A = [1,4]
        Z = [1,2]
        X = [0.75,0.25]

    return X,A,Z



def get_element_ZA(e):

    if (not isinstance(e, basestring)):
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
            "C":, 6,
            "O":, 8
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
