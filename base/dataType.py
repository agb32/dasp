#dasp, the Durham Adaptive optics Simulation Platform.
#Copyright (C) 2004-2016 Alastair Basden and Durham University.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
import types
class dataType:
    """dataType abstract class used to store details about a required variable or object.
    
    Class variables (important to simulation programmer):
     - description - string, description of the type
     - type - Type, the type of variable described
     - dimensions - int or tuple, dimensions of type
     - detail - user defined details
     - val - user defined default value
     - comment - comment for GUI
    @cvar description: A Description of the type
    @type description: String
    @cvar type: The type of the variable
    @type type: e.g. types.IntType, Numeric.ArrayType, "*" etc.
    @cvar dimensions: Dimensions of the object required as a tuple (-1 if don't care)
    @type dimensions: Tuple or int.
    @cvar detail: Extra details describing the object
    @type detail: User defined
    @cvar val: Default value of parameter
    @type val: String
    @cvar comment: Comment for parameter GUI
    @type comment: None or string
    """
    def __init__(self,description=None,typ=None,dimensions=None,detail=None,val=None,comment=None):
        """This class can be used to describe the input/output type required by a module.  If type or
        dimensions don't match, the GUI won't allow the connection.  If
        description doesn't match, the user is warned, but connection is
        allowed.\n
        This is also used by the getParams method to describe the parameters required by the module.
        
        Arguments:
         - description - A text description of the variable.\n
         - typ - The type of the variable, e.g. types.IntType.  Or * if anything will do.\n
         - dimensions - A tuple of dimensions of -1 if don't care (matches any dimensions).\n
         - detail - Any further details, depending on the data type, e.g. Numeric.Float32, None etc.\n

        @param description: Description string
        @type description: String
        @param typ: Object type
        @type typ: Type
        @param dimensions: Dimensions of the object
        @type dimensions: Int or tuple of ints
        @param detail: User defined additional information
        @type detail: User defined
        @param val: Default value of parameter
        @type val: String
        @param comment: comment for parameter GUI
        @type comment: String or None
        """
        self.description=description
        self.type=typ#e.g. types.IntType, Numeric.arrayType, "*" etc
        self.dimensions=dimensions#eg (10,5), -1 (don't care...) etc.
        self.detail=detail#eg None or Numeric.Float32 etc... (any further detail needed)
        self.val=val
        self.comment=comment
    def __repr__(self):
        s="<dataType object %s: %s (%s) %s %s %s>"%(self.description,self.val,self.comment,self.type,self.dimensions,self.detail)
        return s
        
    def compare(self,dataTypes):
        """Compare self with a list of dataTypes.  If dimensions and type
        match with one of these, return it, otherwise return none (no match
        found).

        Arguments:
        dataTypes - a list or (or single) dataTypes objects.
        @param dataTypes: List of (or single) dataTypes objects to be compared with.
        @type dataTypes: List or dataType
        """
        if type(dataTypes)!=types.ListType:
            dataTypes=[dataTypes]
        solution=None
        for dType in dataTypes:
            if dType.type==self.type and dType.dimensions==self.dimensions and self.detail==dType.detail:
                solution=dType
                break
        return solution
