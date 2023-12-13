# Modification for GRDZhadzha

## Makefile
Replace 

```make
GRCHOMBO_SOURCE = ../../Source
```
with
```make
GRDZHADZHA_SOURCE = ../../Source
```

add location of the backgroung simulation parameters

```make
$(GRDZHADZHA_SOURCE)/Background  \
$(GRDZHADZHA_SOURCE)/Matter
```

and copy ```KerrSchildFixedBG.hpp``` from ```Source/FixedBackground``` into ```GRDZhadzha/Source/Background```

remove

```hpp
#include "InitialDataTools.hpp"
```
from ```GRDzhadzha/Source/Background/KerrSchildFixedBG.hpp```

TODO
put ```SpheroidalExtraction.hpp``` in the correct location

why is 

```hpp
emtensor_t
```

reproduced from GRChombo in GRDzhadzha