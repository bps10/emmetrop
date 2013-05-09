#Peripheral photoreceptor activity during accommodation and emmetropization

####Brian P. Schmidt, Maureen Neitz, Jay Neitz

*Neurobiology and Behavior and the Department of Ophthalmology, University of Washington*


### Introduction


####Standard thinking:

* The eye grows until images become sharply focused on the retina
    * Diffusion lenses cause the eye to elongate past the point of emmetropia, resulting in the image coming to focus in front of the eye - myopia.
    * This leads to the logical conclusion that the eye grows until images are sharply focused on the retina.
* The periphery produces the strongest signal used by the eye to figure out when to stop growing.
    * Ablation of the fovea does not interrupt proper emmetropization.
    * Therefore the periphery must be critically important in modulating eye growth.


####Our thinking:

* We agree that visual input on the peripheral retina is key.
* We argue that the relative activity of **photoreceptors in the periphery** are responsible for producing the relevant signal.
* We also argue that **the eye grows until images are blurry** in the periphery, not until they are sharp.  We tested this hypothesis by modeling the activity of a photoreceptor in the periphery during development.


### Overview

Our hypothesis that the eye grows until images become blurry in the periphery has been motivated by two observations:

* Reading an iPod
    * We spend very little time looking at infinity.  Instead we spend much of our time looking at near objects: faces, books, objects we manipulate with our hands, etc.
    * Near work and reading has been strongly correlated with myopia development. 
    * When reading, much of the visual scene is often at a distance.  
    * We frequently focus on near objects (~16in), but the periphery is receiving images that are much more distant. 
    * We must consider this situation when thinking about the activity of photoreceptors.

* Flowers 	
    * Here a relatively hyperopic individual (young) focusing on the flower will be under accommodated, resulting in distant objects that are more sharply focused than they should be.
    * As we approach emmetropia, we focus the flower better, causing the distant objects to become blurrier.
    * This signal will change as a function of age (Anderson et al 2009), potentially modulating the relative activity of photoreceptors in the periphery.

The purpose of this work was to incorporate these ideas into a model that predicts the activity of a cone at 10 degrees eccentricity.

* Schematic - overview of the model.
    1. model natural scenes (well known).
    2. model optical state of the eye (well characterized).
    3. model cone receptive field (well known).


### Model

1. Natural scenes. We used a database of natural images to estimate the power spectrum that we expect to reach the eye (Tkacik et al.). We found the mean spatial frequency, \\(\mathbf{f}\\), spectrum of about 200 images, and fit a power law - a method well described in the scene statistics literature:
	$$s(\mathbf{f}) \, = \, \frac{1}{\mathbf{f}^{\alpha}}$$

	We then corrected the power law for microsaccads. The effects of eye drift on spatial frequency has been extensively studied by the Rucci group and is well described by Brownian motion (2D diffusion equations, transformed into the Fourier domain). \\(D\\) is diffusion constant and \\(\omega\\) is angular frequency:
	$$Q(\mathbf{f}, \omega) \, = \, \frac{\mathbf{f} \, D}{\mathbf{f}^{2} \, \frac{D^{2}}{4} + \omega^{2}}$$

	We incorporate this into the power law with:
	$$S(\mathbf{f}) \, = \, \sum_{\omega =1}^{80} s(\mathbf{f}\)Q(\mathbf{f},\,\omega)$$
	The net effect of eye movements can be seen as a whitening of the spectrum (think white noise, random, equal probability) out to about 10 cycles per degree.

2. Schematic eye. We ray traced the Navarro 1985, 1999 model to derive the modulation transfer functions for various accommodation states.  The Navarro model accounts for accommodation, and off-axis aberrations, but does not include a gradient refractive index (GRIN) lens. This is a limitation, but a minor one as we expect changes to be roughly proportional. *Post-hoc note: the y-axis label of figure 2 should read modulation transfer function*

3. Estimate the spatial frequency content reaching the photoreceptors by combining the functions derived in 1 and 2 above:

	$$r(\mathbf{f}) \, = \, S(\mathbf{f}) \, g(\mathbf{f})$$

4. Cone transfer function, \\(c(\mathbf{f})\\). Fourier transform of the Diff of Gauss receptive field for a cone at 10 degrees eccentricity (2 arcmin spacing between cones).
5. Finally, we can estimate the activity of a cone at 10 degrees eccentricity by multiplying the retinal image, \\(r(\mathbf{f})\\), and the cone transfer function, \\(c(\mathbf{f})\\), and integrating across all spatial frequencies.

$$A \; = \; \int_{0}^{F} r(\mathbf{f}) \, c(\mathbf{f}) \; d\mathbf{f}$$



**Accommodative lag** - changes as a function of age.  We posit this change in lag, which will modulate relative activity of cone photoreceptors in the periphery, could be used by the eye to figure out when to stop growing.

**Reminder**:  Diffraction limit is computed as follows:

$$d \; = \; \frac{\lambda}{2(n \sin{\theta})}$$

where \\(\lambda\\) is the wavelength traveling through the medium, \\(n\\) is the refractive index and \\(\theta\\) is the half angle the cone of light that will be formed at the spot of convergence of the rays - dependent upon the pupil diameter. The denominator is the numerical aperture.


### Results

* Accommodation to near objects results in a significant loss of medium and high spatial frequencies for images of distant objects in the peripheral retina relative to the fovea, reducing the relative activity of photoreceptors there.
* This loss of frequency content is partially ameliorated by accommodative lag that has been observed in young children but decreases during emmetropization. 
* The outcome of this analysis (Table 1) indicates that the accommodation state of the eye substantially modulates the mean activity of a cone photoreceptor.

The summary of the results are shown below.  Object distance and accommodation are represented in Diopters (D) of accommodative demand (computed using the thin lens equation). Cone response is computed by integrating the lines in plot 5 and normalizing by the diffraction limited case (maximum possible under this optical system).

| object | accommodation | cone response|
|:------:|:-------------:|:------------:|
| 0D     | 0D            | 0.89         |
| 2.46D  | 2.46D         | 0.92         |
| 0.16D  | 2.46D         | 0.20         |
| 0.16D  | 0.67D         | 0.43         |

**Conditions**


> 1. Image at infinity, Focus on infinity.
> 2. Image near (16in, book), Focus perfect (16in, book)
> 3. Image far (20ft, background), Focus near (16in, book)
> 4. Image far (20ft, background), Focus underaccommodated (1.79D) (4.9ft)



**Hypotheses**

The table below presents a comparison of the classical idea that the eye grows to clarity with our theory that the eye grows to blurriness.  We present a few important observations that have been made in the myopia literature and the corresponding hypotheses to explain the findings from both theories.

| observation      | eye grows to clarity | eye grows to blurriness |
|:----------------:|:--------------------:|:----------------------:|
| reading causes myopia | ? | high contrast text produces large differences in activity between adjacent cones |
| opsin mutations cause myopia | ? | the mutation increases contrast in cone mosaic |
| emmetropia proceeds normally after fovea ablation | this is contrary to our analysis which shows that the periphery "sees more blur" as the eye grows | the periphery "sees more blur" as the eye grows |
| form deprivation causes myopia| the eye elongates because images are blurred | absence of contrast increases gain so that photoreceptor noise results in large signal differences |


**Notes**

* Opsin mutations cause myopia - Bornholm eye disease.  Cones expressing mutant opsin are rendered inactive, but remain structurally intact.  This leads to a situation where healthy cones expressing a normal opsin sequence sit next to cones that do not transduce light (seen as empty spaces in adaptive optics images), causing a reduction in the lateral inhibition feeding back onto the healthy cone via H2 horizontal cell.  The consequence is a cell producing high activity (contrast), which signals to the eye to continue growing.
* Form deprivation - diffusion lenses decrease the contrast reaching the eye so drastically that the photoreceptors have very low relative activity (they need contrast to be driven). This situation causes the nervous system to crank up the gain - a common principle of early sensory processing.  However, since there is no visually driven signal, the photoreceptors would be amplifying noise - a signal that would artificially look like a high contrast signal.
* Using the thin lens equation to compute accommodative demand: \\(\frac{1}{O} + \frac{1}{I} = \frac{1}{f}\\)
* Change in accommodative lag at 2.46D is not likely to be huge during development, but we do still expect that there are changes and we are only interested in the relative differences.



### Conclusions

* Contrary to the common intuition that during emmetropization the eye grows until images are clear, we compute that considering the statistical environment and the optical transfer functions characteristic of common accommodation states the peripheral retina systematically “sees more blur” and **photoreceptor activity in
the peripheral retina decreases as the eye goes through emmetropization**. 
* Thus, it is fined grained, **sharply focused images which produce high amounts of photoreceptor activity in the periphery** that signal the eye to grow during development and the reduction in that activity that occurs as the eye approaches emmetropia that is responsible for the secession of eye elongation.