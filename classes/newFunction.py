from fenics import *

class NewFunction(Function):

    store = False

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.store = kwargs.get("store")
    
    


        





    


    
