//
//  File.cpp
//  test
//
//  Created by Brian Schmidt on 10/8/12.
//

#include <math.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#include <GL/glut.h>

#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>

#include <Goptical/Math/Vector>
#include <Goptical/Math/VectorPair>

#include <Goptical/Material/Abbe>
#include <Goptical/Material/Base>

#include <Goptical/Data/Set>
#include <Goptical/Data/PlotData>
#include <Goptical/Data/DiscreteSet>
#include <Goptical/Data/Plot>

#include <Goptical/Sys/System>
#include <Goptical/Sys/OpticalSurface>
#include <Goptical/Sys/Source>
#include <Goptical/Sys/SourceRays>
#include <Goptical/Sys/SourcePoint>
#include <Goptical/Sys/Image>
#include <Goptical/Sys/Stop>

#include <Goptical/Curve/Sphere>
#include <Goptical/Curve/Conic>
#include <Goptical/Shape/Disk>

#include <Goptical/Trace/Tracer>
#include <Goptical/Trace/Result>
#include <Goptical/Trace/Distribution>
#include <Goptical/Trace/Sequence>
#include <Goptical/Trace/Params>

#include <Goptical/Light/SpectralLine>

#include <Goptical/Analysis/RayFan>
#include <Goptical/Analysis/Spot>
#include <Goptical/Analysis/Focus>


#include <Goptical/Io/RendererSvg>
#include <Goptical/Io/RendererOpengl>
#include <Goptical/Io/Rgb>
#include <X11/Xlib.h>

// 1. User options: age. GUI.
// 2. Work out GRIN model.
// 3. MTF, PSF.
// 4. Plot in python.
// 5. Off axis.
// 6. Add spectacle lens option!


using namespace Goptical;

class Eye 
{
    
private:    
    float diopters, pupil_size, age, init_diop;
    int print;
    std::string  model;
    
    Sys::System * sys;
    Trace::Tracer * tracer;
    
    Sys::SourceRays  * source_rays;
    Sys::SourcePoint  *  source_point; 
    
    Sys::OpticalSurface * anterior_cornea;
    Sys::OpticalSurface * posterior_cornea;
    Sys::Stop * pupil;
    Sys::OpticalSurface * anterior_lens;
    Sys::OpticalSurface * posterior_lens;
    Sys::Image * image;
    
    Curve::Conic * anterior_cornea_curve;
    Curve::Conic * posterior_cornea_curve;
    Curve::Conic * anterior_lens_curve;
    Curve::Conic * posterior_lens_curve;
    Shape::Disk * ant_cornea_shape;
    Shape::Disk * post_cornea_shape;
    Shape::Disk * lens_shape;
    Shape::Disk * EyeShape;
    Curve::Sphere * EyeCurve;
    
    Analysis::RayFan * fan;
    
public:

    Eye();
    void set_params(int option, std::string mod );
    void set_params(float diop, float pup, std::string mod, float a );
    void SchematicEye( );
    void EyeTracer( );
    void EyePlots( int option );
    void SpotPlot( int option );
    
    float GetCornealThickness(std::string model);
    float GetAnteriorChamber(std::string model, float age, float diopters);
    float GetLensThickness(std::string model, float age, float diopters);
    float GetAxialLength(std::string model, float age, float diopters);
    float GetVitreousLen(std::string model);
    float FindOpticalPower(int opt);
    float Diopters(int option);
}; 

Eye::Eye() {
    
    sys = NULL;
    tracer = NULL;
    source_rays = NULL;
    source_point = NULL;
    anterior_cornea = NULL;
    posterior_cornea = NULL;
    pupil = NULL;
    anterior_lens = NULL;
    posterior_lens = NULL;
    image = NULL;
    
    anterior_cornea_curve = NULL;
    posterior_cornea_curve = NULL;
    anterior_lens_curve = NULL;
    posterior_lens_curve = NULL;
    ant_cornea_shape = NULL;
    post_cornea_shape = NULL;
    lens_shape = NULL;
    EyeShape = NULL;
    EyeCurve = NULL;
    
    fan = NULL;
}

void Eye::set_params(int option, std::string mod = "dubbelman")
{

    std::cout << "enter lens diopters (mm): " << std::endl;
    std::cin >> diopters;
    std::cout << "enter pupil size (mm): " << std::endl;
    std::cin >> pupil_size;
    std::cout << "enter age (years): " << std::endl;
    std::cin >> age;
    
    //age = a;
    model = mod;
    print = 1.0;
    init_diop = diopters;
}

void Eye::set_params(float diop = 0.0, float pup = 5.0, std::string mod = "dubbelman", float a = 20.0)
{

    diopters = diop;
    pupil_size = pup;
    init_diop = diopters;
    age = a;
    model = mod;
    print = 0.0;
    
}


void Eye::SchematicEye()
{
    // create new system
    sys = new Sys::System();
    
    //**********************************************************************
    // set material refraction indexes. Eventually include GRIN for lens:
    //**********************************************************************
    
    ref<Material::AbbeVd>     cornea_refract =      ref<Material::AbbeVd>::create(     1.367,     56.50);
    
    ref<Material::AbbeVd>   extraOcular_refract =   ref<Material::AbbeVd>::create(     1.3374,    49.61);
    
    ref<Material::AbbeVd>   lens_refract =          ref<Material::AbbeVd>::create(     1.42,     48.00);
    
    ref<Material::AbbeVd>   intraOcular_refract =   ref<Material::AbbeVd>::create(     1.336,     50.90);
    
    //**********************************************************************
    // Define system parameters using Dubbelman 2005 or Navarro 1985:
    //**********************************************************************
    std::transform(model.begin(), model.end(), model.begin(), ::tolower);
    
    float corneal_thickness, anterior_chamber, lens_thickness, axial_length, pupil_rad, cornea_ant_k,
    cornea_post_k, cornea_radius_ant, cornea_radius_post, lens_ant_k, lens_post_k, lens_ant_radius,
    lens_post_radius, vitreous_length;
    
    corneal_thickness = GetCornealThickness(model);
    anterior_chamber = GetAnteriorChamber(model, age, diopters);
    lens_thickness = GetLensThickness(model, age, diopters);
    vitreous_length = GetVitreousLen(model);
    axial_length = GetAxialLength(model, age, diopters);
    pupil_rad = pupil_size / 2.0;
        
    if (model == "dubbelman") 
    {
        cornea_ant_k = 0.82;             // Schwarzschild constant (k).
        cornea_post_k = 0.66;            // Schwarzschild constant (k).
        cornea_radius_ant = 7.87;        // radius of curvature (c).
        cornea_radius_post = 6.40;    // radius of curvature (c).
        
        lens_ant_k = -(4.0 - ( 0.5 * diopters ) );                                
        lens_post_k = -3.0;                                                        
        lens_ant_radius =   1.0 / ( 1.0 / (12.7 - 0.058 * age ) + (0.0077 * diopters) ); 
        lens_post_radius = -1.0 / ( 1.0 / (5.9 -  0.013 * age ) + (0.0043 * diopters) ); 
    }
    
    if (model == "navarro")
    { 

        pupil_rad = pupil_size / 2.0;
        
        cornea_ant_k = -0.26;             // Schwarzschild constant (k).
        cornea_post_k = 0;                // Schwarzschild constant (k).
        cornea_radius_ant = 7.72;        // radius of curvature (c).
        cornea_radius_post = 6.50;        // radius of curvature (c).
        
        lens_ant_k = -3.1316;                                
        lens_post_k = -1.0;                                                        
        lens_ant_radius =   10.2; 
        lens_post_radius = -6.0; 
    }
    
    if (print == 1)
    {
        std::cout << "  " << std::endl;
        
        std::cout << "pupil size (mm): "                << pupil_size           << std::endl;
        std::cout << "corneal thickness (mm): "         << corneal_thickness    << std::endl;
        std::cout << "anterior chamber depth (mm): "    << anterior_chamber     << std::endl;
        std::cout << "lens thickness (mm): "            << lens_thickness       << std::endl;
        std::cout << "vitreous length (mm): "           << vitreous_length      << std::endl;
        std::cout << "axial length (mm): "              << axial_length         << std::endl;
        
        std::cout << "cornea: " << std::endl;
        std::cout << "anterior surface, k (mm): "       << cornea_ant_k         << std::endl;
        std::cout << "posterior surface, k (mm): "      << cornea_post_k        << std::endl;
        std::cout << "anterior radius, c (mm): "        << cornea_radius_ant    << std::endl;
        std::cout << "posterior radius, c (mm): "       << cornea_radius_post   << std::endl;
        
        std::cout << "lens: " << std::endl;
        std::cout << "anterior surface, k (mm): "       << lens_ant_k           << std::endl;
        std::cout << "posterior surface, k (mm): "      << lens_post_k          << std::endl;
        std::cout << "anterior radius, c (mm): "        << lens_ant_radius      << std::endl;
        std::cout << "posterior radius, c (mm): "       << lens_post_radius     << std::endl;
    }
    
    //**********************************************************************
    // Cornea: 
    //**********************************************************************
    
    // update cornea parameters
    anterior_cornea_curve = new Curve::Conic( cornea_radius_ant,  cornea_ant_k);  // radius of curvature (c), Schwarzschild constant (k).
    posterior_cornea_curve = new Curve::Conic(cornea_radius_post, cornea_post_k); // radius of curvature (c), Schwarzschild constant (k).
    ant_cornea_shape = new Shape::Disk(5.8);
    post_cornea_shape = new Shape::Disk(4.9); // cornea diameter in mm
    
    anterior_cornea = new Sys::OpticalSurface(Math::Vector3(0, 0, 0),   // position.
                                              *anterior_cornea_curve,   // curve.
                                              *ant_cornea_shape,        // aperture shape.
                                              Material::none,           // material to left.
                                              cornea_refract);          // material to right.
    
    posterior_cornea = new Sys::OpticalSurface(Math::Vector3(0, 0, corneal_thickness), // position.
                                               *posterior_cornea_curve, // curve.
                                               *post_cornea_shape,      // aperature shape.
                                               cornea_refract,          // material to left.
                                               extraOcular_refract);    // material to right.
    
    
    //**********************************************************************
    // add pupil (set location and radius (mm) of pupil (default = 1.5mm).
    //**********************************************************************
    
    pupil = new Sys::Stop(Math::Vector3(0, 0, corneal_thickness + anterior_chamber),  pupil_rad);      
    pupil->set_external_radius( 6.0); // make sure pupil radius is at least as large as cornea.
    
    //**********************************************************************
    // Crystalline lens
    //**********************************************************************
    
    // update lens parameters
    anterior_lens_curve = new Curve::Conic(    lens_ant_radius,     lens_ant_k);     
    posterior_lens_curve = new Curve::Conic(    lens_post_radius,     lens_post_k);             
    lens_shape = new Shape::Disk(6.0); // lens diameter in mm
    
    anterior_lens = new Sys::OpticalSurface(Math::Vector3(0, 0, corneal_thickness + anterior_chamber), 
                                            *anterior_lens_curve,   // curve.
                                            *lens_shape,            // aperture shape.
                                            extraOcular_refract,    // material to left.
                                            lens_refract);          // material to right.
    // position.
    posterior_lens = new Sys::OpticalSurface(Math::Vector3(0, 0, corneal_thickness + anterior_chamber + lens_thickness),     
                                             *posterior_lens_curve, // curve.
                                             *lens_shape,           // aperture shape.
                                             lens_refract,          // material to left.
                                             intraOcular_refract);  // material to right.
    
    //**********************************************************************
    // add all of the optical components.
    //**********************************************************************
    
    sys->add(*anterior_cornea);
    sys->add(*posterior_cornea);
    sys->add(*pupil);
    sys->add(*anterior_lens);
    sys->add(*posterior_lens);
    
    //**********************************************************************
    // set the eye shape:
    //**********************************************************************
    
    EyeCurve = new Curve::Sphere(-10.0); 
    EyeShape = new Shape::Disk(1.0); // eye radius (mm)
    
    image = new Sys::Image(Math::Vector3(0, 0, axial_length), *EyeCurve, *EyeShape); 
    
    sys->add(*image);
    sys->set_entrance_pupil(*anterior_cornea);
    
}


void Eye::EyeTracer()
{
    //**********************************************************************
    // Setup light sources and ray tracer
    //**********************************************************************
    
    source_rays  = new Sys::SourceRays(Math::Vector3(0, 0, -1000));
    source_point = new Sys::SourcePoint(Sys::SourceAtInfinity, Math::vector3_001);  
    
    // add sources to system
    sys->add(*source_rays);
    sys->add(*source_point);
    
    // configure sources
    source_rays->add_chief_rays(*sys);
    source_rays->add_marginal_rays(*sys, pupil_size / 2.0);
    source_rays->add_marginal_rays(*sys, -pupil_size / 2.0);
    
    source_point->clear_spectrum();
    source_point->add_spectral_line(Light::SpectralLine::e);
    
    // ray tracer
    tracer = new Trace::Tracer(*sys);
    
}

void Eye::EyePlots(int option = 1)
{    

            
    switch (option) 
    {
        case 1:
        {
            Io::RendererSvg renderer("eye1.svg", 1200,1200);
            renderer.set_margin_ratio(0.35, 0.25, 0.1, 0.1);
    
            const int number (1);
            // layout plot
            

            renderer.set_page_layout(1,2);
            
            // draw 2d system layout: 
            sys->draw_2d_fit(renderer);
            sys->draw_2d(renderer);
            for (int i = 0; i < number; i++)
            {
                // trace and draw rays from rays source
                sys->enable_single<Sys::Source>(*source_rays);
                tracer->get_trace_result().set_generated_save_state(*source_rays);
                
                renderer.set_page(0);
                tracer->trace();
                tracer->get_trace_result().draw_2d(renderer);
                
                //source_rays.rotate(0, 0.10, 0);
                
                // longitudinal aberration
                fan = new Analysis::RayFan(*sys);
                
                sys->enable_single<Sys::Source>(*source_point);
                ref<Data::Plot> abber_plot = fan->get_plot(Analysis::RayFan::EntranceHeight,
                                                          Analysis::RayFan::LongitudinalDistance);
                
                renderer.set_page(1);
                abber_plot->draw(renderer);
                /*
                int curve_index(0);
                const Data::Set &data = fan->get_plot(Analysis::RayFan::EntranceHeight,
                                                          Analysis::RayFan::LongitudinalDistance)
                                                          ->get_plot_data(curve_index).get_set();
                
                const Data::DiscreteSet * out;
                out = dynamic_cast<const Data::DiscreteSet*> (&data);
                
                std::cout << "out: " << out << std::endl;*/ 
            }
                break;
        }
        case 2:
        {
            Io::RendererSvg renderer("Eye.svg", 1200,800);
            renderer.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer.set_page_layout(1,4);
            for (int i = 0; i < 4; i++)
            {    
                
                set_params(init_diop + (i * 2.0), pupil_size, model, age);
                SchematicEye();
                EyeTracer();
                //EyePlots();
                
                std::cout << "lens (diopters): " << init_diop + (i * 2.0) << "  Age (years): " 
                            << age << "  pupil size (mm):" << pupil_size << std::endl;
                            
                std::cout << "unaccommodated defocus: " << FindOpticalPower(2)
                            << "  best focus (diopter): " << FindOpticalPower(1) << std::endl;
                
                renderer.set_page(i);
                
                sys->enable_single<Sys::Source>(*source_rays);
                tracer->get_trace_result().set_generated_save_state(*source_rays);
                sys->draw_2d_fit(renderer);
                sys->draw_2d(renderer);                
                tracer->trace();
                tracer->get_trace_result().draw_2d(renderer);
                
                
            }
        }    
        
    }
}

void Eye::SpotPlot(int option = 1)
{
    switch (option)
    { 
        case 1:
        {
            const int number (1);
            //  change position of light slightly for a series of plots.
            Io::RendererSvg     renderer("spot.svg",        300 * 1, 300 * number, Io::rgb_black);   
            
            sys->get_tracer_params().set_default_distribution(
                                                Trace::Distribution(Trace::RandomDist, 200)); 
            renderer.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer.set_page_layout(1, number);
            
            for (int i = 0; i < number; i++)
            {
                Analysis::Spot spot(*sys);
                
                renderer.set_page(i);
                spot.draw_diagram(renderer);
                
                //source_point.rotate(0, 0.10, 0);
            }
            break;
        }
        case 2:
        {
            Io::RendererSvg renderer1("spot.svg", 300,1200, Io::rgb_black);
            Io::RendererSvg renderer2("spot_intensity.svg", 640, 480*4);
            Io::RendererSvg renderer3("spotBestFocus.svg", 300,1200, Io::rgb_black);
            Io::RendererSvg renderer4("spot_intensityBestFocus.svg", 640, 480*4);

            renderer1.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer2.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer3.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer4.set_margin_ratio(0.1, 0.1, 0.1, 0.1);

            renderer1.set_page_layout(1,4);
            renderer2.set_page_layout(1,4);
            renderer3.set_page_layout(1,4);
            renderer4.set_page_layout(1,4);

            int init_diop = 0;
            std::ofstream outputfile1, outputfile2;
            outputfile1.open ("Intensity.csv", std::ios::trunc);
            outputfile2.open ("IntensityBestFocus.csv", std::ios::trunc);
            for (int i = 0; i < 4; i++)
            {    
                
                set_params(init_diop + (i * 2.0), pupil_size, model, age);
                SchematicEye();
                EyeTracer();
                EyePlots(1);
                
                renderer1.set_page(i);
                renderer2.set_page(i);
                sys->get_tracer_params().set_default_distribution(
                                            Trace::Distribution(Trace::RandomDist, 500)); 
                
                
                Analysis::Spot spot1(*sys);
                
                spot1.draw_diagram(renderer1);    

                ref<Data::Plot> plot1 = spot1.get_encircled_intensity_plot(100);

                plot1->draw(renderer2);                

                double    radius1 = 0.0; // in mm.
                while ( radius1 < 1.0)
                {
                    outputfile1 << spot1.get_encircled_intensity( radius1 ) << "," << radius1
                        << "," << init_diop + (i * 2.0) <<std::endl;
                    
                    radius1 += 0.001;
                }

        		// now repeat measurements at point of best focus.
                
                renderer3.set_page(i);
                renderer4.set_page(i);
                
                Analysis::Spot spot2(*sys);
                Analysis::Focus focus(*sys);
                image->set_plane(focus.get_best_focus());

                spot2.draw_diagram(renderer3);    

                ref<Data::Plot> plot2 = spot2.get_encircled_intensity_plot(100);

                plot2->draw(renderer4);                

                double    radius2 = 0.0; // in mm.
                while ( radius2 < 1.0)
                {
                    outputfile2 << spot2.get_encircled_intensity( radius2 ) << "," << radius2
                        << "," << init_diop + (i * 2.0) <<std::endl;
                    
                    radius2 += 0.001;
                }

            }
        outputfile1.close();
		outputfile2.close();
        }    
    }
}


float Eye::GetCornealThickness(std::string model)
{    
    float corneal_thickness;
    if (model == "dubbelman")    {corneal_thickness = 0.574;}
    if (model == "navarro")     {corneal_thickness = 0.55;}
    return corneal_thickness; 
}


float Eye::GetAnteriorChamber(std::string model, float age, float diopters)
{
    float anterior_chamber;
    if (model == "dubbelman") // in one paper is written as 3.87 !PLUS! (0.010 ...
    {anterior_chamber = 3.87 - ( 0.010  * age ) - ( diopters * ( 0.048 - 0.0004 * age) );}
    if (model == "navarro") {anterior_chamber = 3.05;}
    
    return anterior_chamber; 
}


float Eye::GetLensThickness(std::string model, float age, float diopters)
{
    float lens_thickness;
    if (model == "dubbelman")
    {lens_thickness = 2.93 + ( 0.0236 * age ) + ( diopters * ( 0.058 - 0.0005 * age));}
    if (model == "navarro") {lens_thickness = 4.0;}
    return lens_thickness;
}


float Eye::GetAxialLength(std::string model, float age, float diopters)
{    
    float a, b, c;
    a = GetAnteriorChamber(model, age, diopters);
    b = GetLensThickness(model, age, diopters);
    c = GetVitreousLen(model);
    return a + b + c;
}


float Eye::GetVitreousLen(std::string model)
{ 
    float vitreous_length;
    if (model == "dubbelman")    {vitreous_length = 16.9935;}
    if (model == "navarro")     {vitreous_length = 16.6;}
    return vitreous_length;
}


float Eye::FindOpticalPower(int opt = 1)
{
    Analysis::Focus        focus(*sys);
    float power, focal_len;
    if (opt == 1) { focal_len = focus.get_best_focus()[0][2]; }
    if (opt == 2) { focal_len = GetAxialLength(model, age, diopters); }
    
    power = 1.0/(focal_len/1000.0);
    return power;
}


float Eye::Diopters(int option = 0)
{
    Analysis::Focus     focus(*sys);
    float relaxed_power, accomm_power, defocus_diopter;
    relaxed_power = FindOpticalPower(2);
    accomm_power = FindOpticalPower(1);
    defocus_diopter = accomm_power - relaxed_power;
    
    if (option == 0)
    {
    std::cout << " " << std::endl;
    std::cout << "focal plane: " << std::endl;
    std::cout << focus.get_best_focus()[0][2] << std::endl;
    std::cout << "defocus diopters: " << std::endl;
    std::cout << defocus_diopter << std::endl;
    return 0.0;
    }
    
    if (option == 1)
    {
    return defocus_diopter;
    }
    
}


int main(int argc, const char * argv[])
{
    int option (0);
    
    for (int i = 1; i < argc; i++) {
        
        if (i == 1) 
        {
            std::string word = argv[1];

            if (word == "plot") {option = 0;}
            if (word == "series") {option = 2;}
            if (word == "loop") {option = 1;} 
                        
        }
        // add in model change param (i.e. Navarro vs Dubbelman) and age option for loop.
        if (i == 2 && option == 0)
        {
            std::string word = argv[2];
            if (argv[i] == "navarro" or argv[i] == "dubbelman") 
                {std::string  mod;
                mod = argv[i];}
            else 
                {std::cout << "sorry model option not understood" << std::endl;}
        }
        //else {std::cout << "model cannot be changed with loop option" << std::endl;}
        
        //else if (i == 2 && option == 1)
        //    std::string word = argv[2];
        //    if (argv[i] == "
        if (i > 2) {std::cout << "sorry additional options not understood" << std::endl;}
    }
    
    if (argc < 1 || option == 0)
    {
    Eye     eye;
    eye.set_params(option);
    eye.SchematicEye();
    eye.EyeTracer();
    eye.EyePlots(1); 
    eye.Diopters();
    }
    
    if (option == 1)
    {

    std::ofstream outputfile;
    outputfile.open ("EyeLSA.csv", std::ios::trunc);
    float AGE, LensAccomm, PupilSize, AccommOptPower, RelaxedOptPower, Defocus_diopters;
    std::string mod;
    
    for (int i = 0; i < 20; i++)
    {
        
        for (int j = 1; j < 21; j++)
        {
            for (int k = 10; k <25; k++)
            {
                LensAccomm = i / 2.0;
                PupilSize = j / 4.0;
                AGE = k;
                
                Eye     eye;
                eye.set_params( LensAccomm, PupilSize, mod = "dubbelman", AGE );
                eye.SchematicEye();
                eye.EyeTracer();
                
                AccommOptPower = eye.FindOpticalPower(1);
                RelaxedOptPower = eye.FindOpticalPower(2);
                Defocus_diopters = eye.Diopters(option);
                
                outputfile << LensAccomm << "," << PupilSize << "," << AGE << "," << AccommOptPower << 
                            "," << RelaxedOptPower << "," << Defocus_diopters << std::endl;
            }
        }
    }
    
    outputfile.close();
     }
    
    if (option ==2)
    {
    std::cout << "starting values... " << std::endl;
    std::cout << "  "  << std::endl;
    
    Eye        eye;
    eye.set_params(option);
    eye.EyePlots(2);
    eye.SpotPlot(2);
    }
    
    
    return 0;
}