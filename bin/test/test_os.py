import os



def stuff():

    return os.path.isfile("stuff")




if stuff():
    print "hi"
else:
    print "na"