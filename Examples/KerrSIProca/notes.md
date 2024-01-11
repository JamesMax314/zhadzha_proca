# Modification for GRDZhadzha

## Makefile
Within the Make file make the following changes:
- Replace ```GRCHOMBO_SOURCE = ../../Source``` with ```GRDZHADZHA_SOURCE = ../../Source```

- add location of the backgroung simulation parameters to 

    ```make
    $(GRDZHADZHA_SOURCE)/Background  \
    $(GRDZHADZHA_SOURCE)/Matter
    ```

## Additional Files

- From GRChombo_Zipeng, copy
    ```bash
    cp GRChombo_Zipeng/Source/FixedBackground/KerrSchildFixedBG.hpp GRDzhadzha/Source/Background
    ```

    ```bash
    cp GRChombo_Zipeng/Source/utils/InitialDataTools.hpp GRDzhadzha/Examples/KerrSIProca
    ```

- SphericalExtraction and SphericalGeometry exist in ```GRChombo/Source/AMRInterpolator``` but not the corresponding spheroidal files:

    ```bash
    cp GRChombo_Zipeng/Source/AMRInterpolator/SpheroidalExtraction.hpp GRDzhadzha/Examples/SpheroidalExtraction
    ```

    ```bash
    cp GRChombo_Zipeng/Source/AMRInterpolator/SpheroidalGeometry.hpp GRDzhadzha/Examples/SpheroidalGeometry
    ```

- There are a number of fixed background files that don't exist in GRDZhadzha:

    ```bash
    cp GRChombo_Zipeng/Source/FixedBackground/FixedBGDensityAndAngularMom.hpp GRDzhadzha/Examples/KerrSIProca/
    ```

    ```bash
    cp GRChombo_Zipeng/Source/FixedBackground/FixedBGEnergyAndAngularMomFlux.hpp GRDzhadzha/Examples/KerrSIProca/
    ```

    ```bash
    cp GRChombo_Zipeng/Source/FixedBackground/FixedBGEvolution.hpp GRDzhadzha/Examples/KerrSIProca/
    ```
- There is a missing utility file:
    ```bash
    cp GRChombo_Zipeng/Source/utils/FluxExtraction.hpp GRDzhadzha/Examples/KerrSIProca/
    ```

- The Coordinate Transformations utility in GRChombo_Zipeng are different from those in GRChombo. Compiling with the GRChombo_Zipeng version fixes this:

    ```bash
    cp GRChombo_Zipeng/Source/utils/CoordinateTransformations.hpp GRDzhadzha/Examples/KerrSIProca/
    ```

## Code Modifications

- In ```GRDzhadzha/Source/Background/ADMFixedBGVars```, 

    ```hpp
    template <class data_t> struct emtensor_t
    ```

    is reproduced from ```GRChombo/Source/CCZ4/CCZ4Geometry.hpp``` and causes a predefined symbol error if not <b>removed</b> in GRDzhadzha.

- The ```GRChombo/Source/AMRInterpolator/SurfaceExtraction.hpp``` file in GRChombo differs in its definition of 

    ```hpp
    struct surface_extraction_params_t
    ```

    vs 

    ```hpp
    template <class SurfaceGeometry> class SurfaceExtraction
    {
    public:
        struct params_t
    }
    ```

    in ```GRChombo_Zipeng/Source/AMRInterpolator/SurfaceExtraction.hpp```

    conesquently, lines 16 and 34 of ```GRDzhadzha/Examples/KerrSIProca/SpheroidalExtraction.hpp``` must be <b>changed</b> from 

    ```hpp
    SurfaceExtraction::params_t
    ```

    to
    
    ```hpp
    surface_extraction_params_t
    ```
