# Command Sequence

**_Assume three layer {0,1,2} data set._**

### Initial Project Folder

>* prms
    * matchparams.txt
    * scriptparams.txt
* mylayout.txt
* topgo.sht

### Create Image Database and 'temp0' Work Folder

```
> ./topgo.sht
> ./dbgo.sht
// ---------------------------------------------
// If using MRC images or if making foldmasks:
>     cd idb    // within idb folder
>    ./fsub.sht 0 2
>    ./freport.sht
>    cd ..      // back to temp0
// ---------------------------------------------
> ./mongo.sht
```

### Extract Same-Layer Point-Pairs

```
> cd temp0
> ./ssub.sht 0 2
> ./sreport.sht 0 2
```

### Solve Montages

```
> ./msub.sht 0 2
> ./mreport.sht 0 2
> ./gathermons.sht 0 2
```

### Create Cross-Layer Work Folder

```
> ./crossgo.sht
```

### Run Coarse Strip-Strip Alignment

```
> cd cross_wkspc
> ./subscapes.sht
> ./lowresgo.sht
```

### Review/Bless Strip-Strip Alignment

* Open LowRes.xml in TrakEM2; fix any layers with `Align/Align Using Manual Landmarks`.
* Save/Close LowRes.xml.

### Create Coarse Scaffold

```
> ./scafgo.sht
```

### Run Block-Block Alignment

```
> ./carvego.sht
> ./bsub.sht 0 2
> ./breport.sht 0 2
```

### Extract Cross-Layer Point-Pairs

```
> cd ..		// now in temp0
> ./dsub.sht 0 2
> ./dreport.sht 0 2
```

### Solve Final Stack

```
> cd stack
> ./runlsq.sht
> ./xviewgo.sht
```

#### Now have aligned TrakEM2 file 'Affine.xml'.

_fin_


