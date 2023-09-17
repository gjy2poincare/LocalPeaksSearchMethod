#                                   LocalPeaksSearchMethod
  **Introduction to the Program: tested on Matlab Version R2023a**
**********************************************************************************
## Program: 
   **Local-Peaks Search Method (LPS Method) - Dispersion equation solving method.**
   **for the following Models:**
*  Free or Fluid-loaded, Single- or Double-Layer, Elastic or Viscoelastic Plates;
*  Free or Fluid-loaded, Single- or Double-Layer, Elastic or Viscoelastic Cylindrical Shells.

## Developers: Jiayuan Gong

## References: 
[1] to be added soon...

## Notes:  
Anyone who optimizes the program, please share wtih the developer.

 ## Model Scheme
* -----------------------------------------------------------------
           Fluid 1 [Fl, Va]: row1,c1
* -----------------------------------------------------------------
       Material 2 [So]: rowvm, Evm0, ytavm, sigmavm      
*  ----------------------------------------------------------------
       Material 1 [So]: rowem, Eem0, ytaem, sigmaem
*  ----------------------------------------------------------------
           Fluid 2 [Fl, Va]: row2, c2
* -----------------------------------------------------------------

## GUI interface Description
**1. Problem Selection**
*  **Model Selection:**     
    (a) PlanarPlate; (b) CylindricalShell
*  **BCs:**    
    Va-free or vacuum; So-solid; Fl-fluid
*  **Mode:**    
    (a) PlanarPlate: S-symmetrical; A-asymmetrical    
    (b) CylindricalShell: L-longitudinal; T-tortional; F-flexural
*  **ncs:**    
    flexural modal order of cylindrical shells
	
**2. Model Size**
*  **Plates:**    
    dem-thickness of material 1;    dvm-thickness of material 2
*  **Cylindrical Shells:**    
    a-inner radius;    b-middle radius;    c-outer radius

**3. Material Parameters**
*  **Fluid 1:**    
    row1, c1 - acoustic parameters of fluid 1
*  **Material 1: elastic:**    
    rowem, Eem0, ytaem, sigmaem ---- parameters of material 1, elastic default
*  **Material 2: viscoelastic:**    
    rowvm, Evm0, ytavm, sigmavm ---- parameters of material 2, viscoelastic default
*  **Fluid 2:**    
    row2, c2 - acoustic parameters of fluid 2

**4. Computation Parameters**
*   **Frequency Range:**     
     f = [fa:df:fb]
*   **Phase Velocity:**    
     cp = [cpa:dcp:cpb]
*   **Attenuation Coefficient:**    
     ki = [kia:dki:kib]
*   **Control Parameters:**    
    (a) kur - used to determine a local-peak    
    (b) err - roots' precision
* **Notes:** wavenumber k = w/cp+i*ki

**5. Inialize & Save**
*  **Initialize:**    
     initialize the parameters using 'Initializer.mat' file saved last time
*  **Save:**    
     save the parameters set currently to 'Initializer.mat' file
	 
**6. Run**
  run the program
  
## Additional Notes
  When the program is running, some dialog boxes will pop-up, which is intended to 
  check the roots obtained through observation by yourself. 
*  **Choice:** to check roots manually?    
      please click 'Yes' or 'No' to determine whether check roots yourself.
*  **Check roots:** Begin to check roots?     
      'Yes'-check; 'No'-not check; 'Never'-not check any more
*  **Is this a root?**     
     through observation to determine a root or not, press "ctrl+C" to changer the current figure    
          'Yes'-is a root; 'No'-not a root; 'Cancel? not sure

