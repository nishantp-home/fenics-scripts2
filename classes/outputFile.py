from fenics import *

class OutputFile:
    def __init__(self, fileName, listOfFields) -> None:
        self._fileName = fileName
        self.listOfFields = listOfFields
        self._fileObject = self._openFile()

    def _openFile(self):
        fileNameWithExtension = self._fileName + ".xdmf"
        fileObject = XDMFFile(fileNameWithExtension)  
        fileObject.parameters.update({"functions_share_mesh": True,"rewrite_function_mesh": False})
        return fileObject
        
    def exportFile(self, time=0):
        for field in self.listOfFields.values():
            if field.store == True:
                self._fileObject.write(field, time)






        