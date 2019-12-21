import os

class SedonaParam(dict):

    def __init__(self,template_file=None,defaults_file=None):

        self.stringlist = []
        self.helper_lines = ""
        self.param_dict = {}

        if (defaults_file is None):
            defaults_file = (os.environ['SEDONA_HOME']) + "/defaults/sedona_defaults.lua"

        self.defaults_file = defaults_file
        self.read_defaults_file(self.defaults_file)
        self.defaults_dict = self.param_dict.copy()
        self.defaults_file = "\"" + self.defaults_file + "\""

        if (template_file is not None):
            self.set_template(template_file)

    #    self.param_dict.update({"defaults_file":self.defaults_file})



    def __getitem__(self, key):
        return self.param_dict[key]


    def __setitem__(self, key, item):

        if (key not in self.defaults_dict.keys()):
            import difflib
            message = key + " is not a valid Sedona parameter\n"
            message += "Did you mean one of: \n"
            for e in difflib.get_close_matches(key, self.defaults_dict.keys()):
                message += "> " + str(e) + "\n"
            raise ValueError(message)

        if (key in self.stringlist):
            self.param_dict[key] = "\"" + item + "\""
        else:
            self.param_dict[key] = item


    def read_defaults_file(self,fname):

        try:
            fin = open(fname,"r")
        except ValueError:
            print ("Can't open defaults file " + fname + "\n")

        for line in fin:

            # remove comments, tabs, etc...
            line_clean = line.partition("--")[0]
            line_clean = line_clean.translate(None, '\t\n ')

            if (not "=" in line_clean):
                continue

            data = line_clean.split("=")
            if (len(data) != 2):
                continue
            self.param_dict.update({data[0]:data[1]})
            if ("\"" in data[1]):
                self.stringlist.append(data[0])

        fin.close()




    def set_template(self,fname=None):


        tm_dir = os.environ['SEDONA_HOME'] + "/defaults/templates"

        if (fname is None):
            print("\nAvailable parameter templates are: ")

            templates = []

            # get all the template names
            for subdir, dirs, files in os.walk(tm_dir):
                for file in files:
                    templates.append(os.path.join(subdir, file))

#            for filename in os.listdir(tm_dir):
#                if ("template_" not in filename):
#                    continue
#                templates.append(filename)

            for i,t in enumerate(templates):
                name = t.replace(tm_dir + "/","")
                name = name.replace(".lua","")
                print str(i) + ') ' + name

            num = raw_input("choose template number >")

            num = int(num)
            if (num < 0 or num >= len(templates)):
                print "Not a valid selection"
                return
            fname = templates[num]
        else:
            fname = tm_dir + "/" + fname + ".lua"


        try:
            fin = open(fname,"r")
        except ValueError:
            print ("Can't open template file " + tm_dir + templates[num] + "\n")

        for line in fin:

            # remove comments, tabs, etc...
            line_clean = line.partition("--")[0]
            line_clean = line_clean.translate(None, '\t\n ')

            if (not "=" in line_clean):
                continue

            data = line_clean.split("=")
            if (len(data) != 2):
                continue

            if (data[0] != "defaults_file"):
                self.param_dict[data[0]] = data[1]
            else:
                self.defaults_file = data[1]
        fin.close()

    def __str__(self):


        line = "----------------------------------\n"
        lasthead = ""


        for key in sorted(self.param_dict.keys()):
            if (key not in self.defaults_dict.keys()):
                line += "{0:35} = {1}\n".format(key,self.param_dict[key])
        line += "----------------------------------\n"

        for key in sorted(self.param_dict.keys()):
            if (key in self.defaults_dict.keys()):
                if (self.param_dict[key] != self.defaults_dict[key] or all==True):
                    head = (key.split("_"))[0]
                    if (head != lasthead):
                        line += "--\n"
                        lasthead = head

                    line += "{0:35} = {1}\n".format(key,self.param_dict[key])

        line += "----------------------------------\n"
        return line

    def write(self,filename="param.lua",all=False):


        for key in sorted(self.param_dict.keys()):
            if ("\"" in key):
                print "string",key,self.param_dict[key]

        try:
            fout = open(filename,"w")
        except:
            print "Can't open output param filename " + filename

        linebreak = "--------------------------------------------\n"

        # print user define variables
        line = linebreak
        for key in sorted(self.param_dict.keys()):
            if (key not in self.defaults_dict.keys()):
                line += "{0:35} = {1}\n".format(key,self.param_dict[key])
        line += linebreak
        fout.write(line)

        # print defaults file
        line = "{0:35} = {1}\n".format("defaults_file",self.defaults_file)
        line += linebreak
        fout.write(line)


        # print out a few params to the top
        header_keys = ["model_file","grid_type"]
        line = ""
        for key in header_keys:
            if (self.param_dict[key] != self.defaults_dict[key] or all==True):
                line += "{0:35} = {1}\n".format(key,self.param_dict[key])
        line += linebreak
        fout.write(line)

        lasthead = ""
        for key in sorted(self.param_dict.keys()):
            if (key in header_keys):
                continue
            if (key in self.defaults_dict.keys()):
                if (self.param_dict[key] != self.defaults_dict[key] or all==True):
                    head = (key.split("_"))[0]
                    if (head != lasthead):
                        fout.write("--\n")
                        lasthead = head

                    line = "{0:35} = {1}\n".format(key,self.param_dict[key])
                    fout.write(line)
        fout.write(linebreak)
