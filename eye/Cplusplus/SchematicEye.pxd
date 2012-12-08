# distutils: language = c++
# distutils: sources = SchematicEye.cc

# distutils: include_dirs = c:\MinGW\usr\local\include
# distutils: libraries = Goptical

cdef extern from "SchematicEye.h":

    cdef cppclass Eye:
        Eye() except +

        
        void set_params(float diop, float pup, str mod, float a )
        void SchematicEye( )
        void EyeTracer( )
        void EyePlots( int option )
        void SpotPlot( int option )
        
        float GetCornealThickness(str model)
        float GetAnteriorChamber(str model, float age, float diopters)
        float GetLensThickness(str model, float age, float diopters)
        float GetAxialLength(str model, float age, float diopters)
        float GetVitreousLen(str model)
        float FindOpticalPower(int opt)
        float Diopters(int option)

