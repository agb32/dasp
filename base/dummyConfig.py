class config:
    """A dummy config, that can be used for testing and debugging"""
    def __init__(self,values):
        """values is a dictionary of values, with keys being the name of the variable to be searched for"""
        self.val=values#a dictionary of values
    def getVal(self,name,default="nodefault"):
        if self.val.has_key(name):
            val=self.val[name]
        elif default!="nodefault":
            val=default
        else:
            raise Exception("Variable %s not found in dummy config"%name)
        return val
    def setSearchOrder(self,val1=None,val2=None):
        pass
