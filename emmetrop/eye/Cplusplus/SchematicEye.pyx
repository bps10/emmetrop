
cdef class PyEye:
    cdef Eye* thisptr # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisprt = new Eye( )
    def __delaloc__(self):
        del self.thisptr
    def set_params(self, float diop, float pup, str mod, float a ):
        self.thisptr.set_params(diop, pup, mod, a)
    def SchematicEye(self):
        self.thisptr.SchematicEye()  
    def EyeTracer(self):
        self.thisptr.EyeTracer()
    def EyePlots(self, int option):
        self.thisptr.EyePlots(option)
    def SpotPlot(self, int option):
        self.thisptr.SpotPlot(option)
    def GetCornealThickness(self, str model):
        self.thisptr.GetCornealThickness(model)
    def GetAnteriorChamber(self, str model, float age, float diopters):
        self.thisptr.GetAnteriorChamber(model, age, diopters)
    def GetLensThickness(self, str model, float age, float diopters):
        self.thisptr.GetLensThickness(model, age, diopters)
    def GetAxialLength(self, str model, float age, float diopters):
        self.thisptr.GetAxialLength(model, age, diopters)
    def GetVitreousLen(self, str model):
        self.thisptr.GetVitreousLen(model)
    def FindOpticalPower(self, int opt):
        self.thisptr.FindOpticalPower(opt)
    def Diopters(self, int option):
        self.thisptr.Diopters(option)
        
print PyEye   