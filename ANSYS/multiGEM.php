<?php

$VOLUMES = array();
$ERROR = ''; // Error string

function println($str) {
    print($str.PHP_EOL);
}

function printcmt($str) {
    print(PHP_EOL.'! '.$str.PHP_EOL);
}

function setV($name) {
    global $VOLUMES;
    $loop=true;
    $i=1;
    while($loop) {
        
        if(!array_key_exists($i, $VOLUMES)) {
            $VOLUMES[$i] = $name;
            break;
        }
        $i++;
    }
    return $i;
}

function delV($id) {
    global $VOLUMES, $ERROR;
    if(array_key_exists($id, $VOLUMES)) unset($VOLUMES[$id]);
    else $ERROR .= 'delV ID='.$id.' not found<br />';
}

function generateFile($GEM, $gemLayers) {
    
	global $VOLUMES, $ERROR;

	// Total amount of pieces of one GEM hole
	$gemElements = 5;                   

	// Total thickness of one GEM layer
	$layerThickness = 2*$GEM['metal'] + $GEM['kapton'];

	// Total detector thickness
	$totalThickness = $GEM['DRIFTT'] + $GEM['INDUCTT'] + $gemLayers*$layerThickness;
	for($i=1; $i<$gemLayers; $i++) {
		$totalThickness = $totalThickness + $GEM['TRANSFERT'.$i];
	}

	// BEGIN FILE GENERATION

    echo '<pre>';
    
    ob_start();
    
    println('FINISH');
    println('/CLEAR,START');
    println('/PREP7');  

    printcmt('No polynomial elements');
    println('/PMETH,OFF,1');

    printcmt('Set electric preferences');
    println('KEYW,PR_ELMAG,1');
    println('KEYW,MAGELC,1');   

    printcmt('Select element');
    println('ET,1,SOLID123');

    printcmt('Construct the GEM');
    println('pitch = '.$GEM['pitch']);
    println('kapton = '.$GEM['kapton']);
    println('metal = '.$GEM['metal']);
    println('outdia = '.$GEM['outdia']);
    println('middia = '.$GEM['middia']);
    println('rim  = '.$GEM['rim']);
    println('drift = '.$GEM['DRIFTT']);
    for($i=1; $i<$gemLayers; $i++) {
        println('transfer'.$i.' = '.$GEM['TRANSFERT'.$i]);
    }
    println('induct = '.$GEM['INDUCTT']);

    /*
     * Materials
     */

    printcmt("Material properties");
    println("MP,PERX,1,1e10");  // Metal permitivity
    println("MP,RSVX,1,0.0");   // Metal resisitivty
    println("MP,PERX,2,1.0");   // Gas
    println("MP,PERX,3,4.0");   // Permittivity of kapton

    /*
     * PART 1: make the blocks
     */


    // GEM LAYERS
    $zIndex = $GEM['DRIFTT'];
    for($i=0; $i<$gemLayers; $i++) {

        printcmt('GEM LAYER '.($i+1));

        // LOWER METAL
        println("BLOCK, 0, ".$GEM['pitch']."/2, 0, sqrt(3)*".$GEM['pitch']."/2, ".$zIndex.", ".($zIndex+$GEM['metal']));
        setV("METAL LOWER LAYER=".($i+1));
        $zIndex+=$GEM['metal'];

        // KAPTON
        println("BLOCK, 0, ".$GEM['pitch']."/2, 0, sqrt(3)*".$GEM['pitch']."/2, ".$zIndex.", ".($zIndex+$GEM['kapton']));
        setV("KAPTON LAYER=".($i+1));
        $zIndex+=$GEM['kapton'];

        // UPPER METAL
        println("BLOCK, 0, ".$GEM['pitch']."/2, 0, sqrt(3)*".$GEM['pitch']."/2, ".$zIndex.", ".($zIndex+$GEM['metal']));
        setV("METAL UPPER LAYER=".($i+1));
        $zIndex+=$GEM['metal'];

        // TRANSFER GAP (but not for the last layer!)
        if($i != $gemLayers-1) $zIndex+=$GEM['TRANSFERT'.($i+1)];

    }
    $zIndex += $GEM['INDUCTT'];

    // GAS GAP
    printcmt('TOTAL GAS GAP');
    println("BLOCK, 0, ".$GEM['pitch']."/2, 0, sqrt(3)*".$GEM['pitch']."/2, 0, ".$totalThickness);
    setV("GAS_GAP");

    /*
     * PART 2: create the cones and cylinders in each GEM layer
     */

    // Make cones and cylinders with x=y=0 LABEL=(1)
    for($i=0; $i<$gemLayers; $i++) {

        printcmt("Make cut-out pieces in layer ".($i+1)." at x=y=0");

        // move $zIndex to the middle of the layer
        $zIndex = $GEM['DRIFTT'] + $GEM['metal'] + 0.5*$GEM['kapton'] + $i*$layerThickness;
        for($j=1; $j<($i+1); $j++) {
            $zIndex = $zIndex + $GEM['TRANSFERT'.$j];
        }
        println("WPOFFS, 0, 0, ".$zIndex);

        println("CONE, outdia/2, middia/2, -kapton/2,   0, 0, 360");
        $v1=setV('CONE (1) LOWER LAYER='.($i+1));

        println("CONE, middia/2, outdia/2, 0, kapton/2, 0, 0, 360");
        $v2=setV('CONE (1) UPPER LAYER='.($i+1));

        println("WPOFFS, 0, 0, kapton/2");
        println("CYL4,   0, 0, rim/2, ,,, metal");
        $v3=setV('CYLINDER (1) UPPER LAYER='.($i+1));

        println("WPOFFS, 0, 0, -kapton");
        println("CYL4,   0, 0, rim/2, ,,, -metal");
        $v4=setV('CYLINDER (1) LOWER LAYER='.($i+1));
        println("WPOFFS, 0, 0, kapton/2"); // return to center of layer

        // Reset offset to origin
        println("WPOFFS, 0, 0, -".$zIndex);

        // Now, we have created 4 volumes with numbers $v1, $v2, $v3 and $v4
        // We merge these volumes into one a new volume
        // The others are deleted!
        println("VADD, ".$v1.", ".$v2.", ".$v3.", ".$v4."");
        setV("HOLE (1) LAYER=".($i+1));
        delV($v1);
        delV($v2);
        delV($v3);
        delV($v4);
    }

    // Make cones and cylinders with x=pitch/2, y=sqrt(3)*pitch/2 LABEL=(2)
    for($i=0; $i<$gemLayers; $i++) {

        printcmt("Make cut-out pieces in layer ".($i+1)." at x=pitch/2, y=sqrt(3)/2*pitch");

        // move $zIndex to the middle of the layer
        $zIndex = $GEM['DRIFTT'] + $GEM['metal'] + 0.5*$GEM['kapton'] + $i*$layerThickness;
        for($j=1; $j<($i+1); $j++) {
            $zIndex = $zIndex + $GEM['TRANSFERT'.$j];
        }
        println("WPOFFS, pitch/2, sqrt(3)*pitch/2, ".$zIndex);

        println("CONE, outdia/2, middia/2, -kapton/2, 0, 0, 360");
        $v1=setV('CONE (2) LOWER LAYER='.($i+1));

        println("CONE, middia/2, outdia/2, 0, kapton/2, 0, 0, 360");
        $v2=setV('CONE (2) UPPER LAYER='.($i+1));

        println("WPOFFS, 0, 0, kapton/2");
        println("CYL4,   0, 0, rim/2, ,,, metal");
        $v3=setV('CYLINDER (2) UPPER LAYER='.($i+1));

        println("WPOFFS, 0, 0, -kapton");
        println("CYL4,   0, 0, rim/2, ,,, -metal");
        $v4=setV('CYLINDER (2) LOWER LAYER='.($i+1));
        println("WPOFFS, 0, 0, kapton/2"); // return to center of layer

        // Reset offset to origin
        println("WPOFFS, -pitch/2, -sqrt(3)*pitch/2, -".$zIndex);    

        // Now, we have created 4 volumes with numbers $v1, $v2, $v3 and $v4
        // We merge these volumes into one a new volume
        // The others are deleted!
        println("VADD, ".$v1.", ".$v2.", ".$v3.", ".$v4."");
        setV("HOLE (2) LAYER=".($i+1));
        delV($v1);
        delV($v2);
        delV($v3);
        delV($v4); 
    }


    /*
     * PART 3: substract the cones and cylinders from the blocks
     */

    $totalLayers = 1 + $gemLayers*(1 + 1 + 1); // gas gap, metal, kapton, metal
    $totalHoles = $gemLayers*2;

    for($i=0; $i<$gemLayers; $i++) {

        $layerIndex = 1 + $i*(1 + 1 + 1);
        $holeIndex = $totalLayers + 5;

        printcmt("Substract holes in layer ".($i+1));

        // Lower metal
        println("VSBV,  ".$layerIndex.",  ".($holeIndex+$i).", , DELETE, KEEP"); // Becomes volume TMP, frees
        $tmp = setV("TMP");
        println("VSBV,  ".$tmp.",  ".($holeIndex+$totalHoles/2+$i).", , DELETE, KEEP");
        delV($tmp);

        // Kapton
        println("VSBV,  ".($layerIndex+1).",  ".($holeIndex+$i).", , DELETE, KEEP"); // Becomes volume TMP, frees
        $tmp = setV("TMP");
        println("VSBV,  ".$tmp.",  ".($holeIndex+$totalHoles/2+$i).", , DELETE, KEEP");
        delV($tmp);

        // Upper metal (and delete the holes)
        println("VSBV,  ".($layerIndex+2).",  ".($holeIndex+$i).", , DELETE, DELETE"); // Becomes volume TMP, frees
        $tmp = setV("TMP");
        println("VSBV,  ".$tmp.",  ".($holeIndex+$totalHoles/2+$i).", , DELETE, DELETE");
        delV($tmp);

        // Delete volumes
        delV($holeIndex+$i);
        delV($holeIndex+$totalHoles/2+$i);
    }


    printcmt("Subtract the kapton and metal from the gas");
    $tmp = count($VOLUMES); // This is the current GAS ID
    for($i=0; $i<$gemLayers; $i++) {

        $layerIndex = 1 + $i*(1 + 1 + 1);

        // Gas ID = 1, becomes $tmp
        println("VSBV, ".$tmp.", ".$layerIndex.", , KEEP, KEEP");
        $tmp1 = setV("TMP");
        //delV($tmp);
        println("VSBV, ".$tmp1.", ".($layerIndex+1).", , KEEP, KEEP"); 
        $tmp2 = setV("TMP");
        //delV($tmp1);
        println("VSBV, ".$tmp2.", ".($layerIndex+2).", , KEEP, KEEP");
        $tmp3 = setV("GAS_GAP"); // This is the new final gas gap
        //delV($tmp2);

        // delete the empty volume IDs
        // The result is that $tmp3 becomes $tmp, i.e. the original setup!
        println("VDEL, ".$tmp);
        println("VDEL, ".$tmp1);
        println("VDEL, ".$tmp2);
        println("NUMCMP, VOLU"); 
        delV($tmp);
        delV($tmp1);
        delV($tmp2);
        delV($tmp3);
        setV("GAS_GAP");
    }

    /*
     * PART 4: glue, paint the different pieces and assign the materials
     * We only have to glue the metal-kapton layers, because the metal-gas and kapton-gas are already glued due to the substraction! (see manual..)
     */

    // Gluing metal-kapton-metal for each layer
    // Because we glue for each layer two elements, the volume ID ordering remains the same
    // Thus, we do not make any changes to the $VOLUME array
    printcmt('Gluing the pieces metal-kapton-metal');
    for($i=0; $i<$gemLayers; $i++) {

        $layerIndex = 1 + $i*(1 + 1 + 1);

        println("VGLUE, ".($layerIndex).", ".($layerIndex+1)); // It deletes the highest number
        $tmp = setV("TMP"); // new volume is stored in $tmp
        println("VGLUE, ".$tmp.", ".($layerIndex+2)); 
        delV($tmp);
    }

    printcmt('Paint the pieces');
    println("/COLOR, VOLU, RED, ".  count($VOLUMES)); // gas
    for($i=0; $i<$gemLayers; $i++) {

        $layerIndex = 1 + $i*(1 + 1 + 1);

        println("/COLOR, VOLU, BLACK, ".($layerIndex));       // Lower metal
        println("/COLOR, VOLU, ORANGE, ".($layerIndex+1));    // Kapton
        println("/COLOR, VOLU, BLACK, ".($layerIndex+2));     // Upper metal
    }


    printcmt('Assign materials');
    println("VSEL, S, , , ".count($VOLUMES)); // gas
    println("VATT, 2, ,1");
    for($i=0; $i<$gemLayers; $i++) {

        $layerIndex = 1 + $i*(1 + 1 + 1);

        println("VSEL, S, , , ".($layerIndex));       // Lower metal
        println("VATT, 1, ,1");
        println("VSEL, S, , , ".($layerIndex+1));     // Kapton
        println("VATT, 3, ,1");
        println("VSEL, S, , , ".($layerIndex+2));     // Upper metal
        println("VATT, 1, ,1");
    }


    /*
    // Gluing metal-kapton-metal with gas 
    $t=array();
    for($i=0; $i<$gemLayers; $i++) {

        $layerIndex = 1 + 1 + $i*(1 + 1 + 1);

        println("VGLUE, 1, ".($layerIndex)); // the lowest number (i.e. 1, the gas) is kept, it deletes the highest number (i.e. $layerIndex)
        delV($layerIndex);
        println("BLOCK, 0, 0.1, 0, 0.1, 0,-0.1");
        array_push($t, setV("TEST".$layerIndex)); // occupy 
        setV("TMP".$layerIndex);

        println("VGLUE, 1, ".($layerIndex+1));
        delV(($layerIndex+1));
        println("BLOCK, 0, 0.1, 0, 0.1, 0,-0.1");
        array_push($t, setV("TEST".($layerIndex+1))); // occupy
        setV("TMP".($layerIndex+1));

        println("VGLUE, 1, ".($layerIndex+2));
        delV(($layerIndex+2));
        println("BLOCK, 0, 0.1, 0, 0.1, 0,-0.1");
        array_push($t, setV("TEST".($layerIndex+2))); // occupy
        setV("TMP".($layerIndex+2));

    }

    foreach ($t as $value) {

        println("VDELE, ".$value);
    }
    println("NUMCMP, VOLU");
    */

    printcmt("Voltage boundaries on the drift and induction plane");
    println("ASEL, S, LOC, Z, 0");
    println("DA, ALL, VOLT, 0"); // zero voltage on first plane
    println("ASEL, S, LOC, Z, ".$totalThickness);
    $v = $GEM['DRIFTV'] + $GEM['INDUCTV'];
    for($i=1; $i<$gemLayers; $i++) {
        $v = $v + $GEM['TRANSFERV'.$i] + $GEM['VGEM'.$i];
    } 
    $v = $v + $GEM['VGEM'.$gemLayers];
    println("DA, ALL, VOLT, ".$v);
    
    // Symmetric boundary conditions GAS
    printcmt("Symmetric boundary conditions on the sides: GAS");
    end($VOLUMES);
    $gasID = key($VOLUMES);
    println("VSEL, S, , , ".$gasID); // Select gas volume
    println("ASLV, S"); // Select all the areas from the selected volumes
    println("ASEL, R, LOC, X, 0"); // Select the area at X=0
    println("DA, ALL, SYMM"); // Apply DOF constraints on selected areas (ALL)
    
    println("VSEL, S, , , ".$gasID); // Select gas volume again.. repeat all the steps
    println("ASLV, S");
    println("ASEL, R, LOC, X, ".$GEM['pitch']."/2");
    println("DA, ALL, SYMM");
    
    println("VSEL, S, , , ".$gasID); 
    println("ASLV, S");
    println("ASEL, R, LOC, Y, 0");
    println("DA, ALL, SYMM");
    
    println("VSEL, S, , , ".$gasID);
    println("ASLV, S");
    println("ASEL, R, LOC, Y, sqrt(3)*".$GEM['pitch']."/2");
    println("DA, ALL, SYMM");
    
    // Symmetric boundary conditions KAPTON
    for($i=0; $i<$gemLayers; $i++) {

        $layerIndex = 1 + $i*(1 + 1 + 1);
        $kaptonID = $layerIndex+1;
        printcmt('Symmetric boundary conditions on the sides: KAPTON (layer: '.($i+1).')');
        
        println("VSEL, S, , , ".$kaptonID);
        println("ASLV, S");
        println("ASEL, R, LOC, X, 0");
        println("DA, ALL, SYMM");

        println("VSEL, S, , , ".$kaptonID);
        println("ASLV, S");
        println("ASEL, R, LOC, X, ".$GEM['pitch']."/2");
        println("DA, ALL, SYMM");

        println("VSEL, S, , , ".$kaptonID); 
        println("ASLV, S");
        println("ASEL, R, LOC, Y, 0");
        println("DA, ALL, SYMM");

        println("VSEL, S, , , ".$kaptonID);
        println("ASLV, S");
        println("ASEL, R, LOC, Y, sqrt(3)*".$GEM['pitch']."/2");
        println("DA, ALL, SYMM");
    } 
    
    // Voltage boundary conditions
    $v = $GEM['DRIFTV'];
    for($i=0; $i<$gemLayers; $i++) {
        
        $layerIndex = 1 + $i*(1 + 1 + 1);      
        if($i != 0) $v = $v + $GEM['TRANSFERV'.$i];
        
        printcmt('Voltage boundary condition on the metal (layer: '.($i+1).')');
        println("VSEL, S, , , ".$layerIndex);
        println("ASLV, S");
        println("DA, ALL, VOLT, ".$v);
        println("VSEL, S, , , ".($layerIndex+2));
        println("ASLV, S");
        $v = $v + $GEM['VGEM'.($i+1)];
        println("DA, ALL, VOLT, ".$v);
    } 
   
    printcmt("Meshing options");
    println("VSEL, S,,, 1, ".count($VOLUMES)); // Select all volumes
    println("ASLV, S"); // Select all areas containing volumes
    println("MSHKEY,0"); // use free meshing
    println("SMRT, 1");

    // Select all kaptons and GAS for meshing
    // The metals are field = 0, so no meshing needed!
    println("VSEL, S, , , ".$gasID);
    for($i=0; $i<$gemLayers; $i++) {

        $layerIndex = 1 + $i*(1 + 1 + 1);
        println("VSEL, A, , , ".($layerIndex+1  ));
    } 

    println("VMESH, ALL");
    
    printcmt("Solve the field");
    println("/SOLU");
    println("SOLVE");
    println("FINISH");

    printcmt("Display the solution");
    println("/POST1"); // Enters the database results postprocessor.
    println("/EFACET, 1"); // Specifies the number of facets per element edge for PowerGraphics displays.
    println("PLNSOL, VOLT, , 0"); // Displays results as continuous contours.

    printcmt("Write the solution to files");
    println("/OUTPUT, PRNSOL, lis");
    println("PRNSOL");
    println("/OUTPUT");
    println("/OUTPUT, NLIST, lis");
    println("NLIST,,,,COORD");
    println("/OUTPUT");
    println("/OUTPUT, ELIST, lis");
    println("ELIST");
    println("/OUTPUT");
    println("/OUTPUT, MPLIST, lis");
    println("MPLIST");
    println("/OUTPUT");
    
    
    
    // WEIGHTNING FIELD
    printcmt("Calculate weighting field (readout)");
    println("/SOLU");
    println("LSCLEAR, ALL");
   
    printcmt("Voltage boundaries on the drift (=0 V) and induction (=1 V) plane");
    println("ASEL, S, LOC, Z, 0");
    println("DA, ALL, VOLT, 0");
    println("ASEL, S, LOC, Z, ".$totalThickness);
    println("DA, ALL, VOLT, 1");

    // Voltage boundary conditions
    for($i=0; $i<$gemLayers; $i++) {

        $layerIndex = 1 + $i*(1 + 1 + 1);

        printcmt('Voltage boundary condition on the metal (layer: '.($i+1).')');
        println("VSEL, S, , , ".$layerIndex);
        println("ASLV, S");
        println("DA, ALL, VOLT, 0");
        println("VSEL, S, , , ".($layerIndex+2));
        println("ASLV, S");
        println("DA, ALL, VOLT, 0");
    }

    printcmt("Meshing options");
    println("VSEL, S,,, 1, ".count($VOLUMES)); // Select all volumes
    println("ASLV, S"); // Select all areas containing volumes
   
    printcmt("Solve the field");
    println("SOLVE");

    printcmt("Display the solution");
    println("/POST1"); // Enters the database results postprocessor.
    println("/EFACET, 1"); // Specifies the number of facets per element edge for PowerGraphics displays.
    println("PLNSOL, VOLT, , 0"); // Displays results as continuous contours.

    printcmt("Write the solution");
    println("/OUTPUT, WSOL, lis");
    println("PRNSOL");
    println("/OUTPUT");


    printcmt("Calculate weighting field (drift)");
    println("/SOLU");
    println("LSCLEAR, ALL");
   
    printcmt("Voltage boundaries on the drift (=1 V) and induction (=0 V) plane");
    println("ASEL, S, LOC, Z, 0");
    println("DA, ALL, VOLT, 1");
    println("ASEL, S, LOC, Z, ".$totalThickness);
    println("DA, ALL, VOLT, 0");

    // Voltage boundary conditions
    for($i=0; $i<$gemLayers; $i++) {

        $layerIndex = 1 + $i*(1 + 1 + 1);

        printcmt('Voltage boundary condition on the metal (layer: '.($i+1).')');
        println("VSEL, S, , , ".$layerIndex);
        println("ASLV, S");
        println("DA, ALL, VOLT, 0");
        println("VSEL, S, , , ".($layerIndex+2));
        println("ASLV, S");
        println("DA, ALL, VOLT, 0");
    }

    printcmt("Meshing options");
    println("VSEL, S,,, 1, ".count($VOLUMES)); // Select all volumes
    println("ASLV, S"); // Select all areas containing volumes
   
    printcmt("Solve the field");
    println("SOLVE");

    printcmt("Display the solution");
    println("/POST1"); // Enters the database results postprocessor.
    println("/EFACET, 1"); // Specifies the number of facets per element edge for PowerGraphics displays.
    println("PLNSOL, VOLT, , 0"); // Displays results as continuous contours.

    printcmt("Write the solution");
    println("/OUTPUT, WSOLD, lis");
    println("PRNSOL");
    println("/OUTPUT");

    $content = ob_get_clean();
    $fh = fopen("output/GEM.inp", 'w') or die("can't open file");
    fwrite($fh, $content);
    fclose($fh);
    echo $content;

    echo '</pre>';

    echo '<h3>Assigned volumes</h3>';
    echo '<pre>';
    print_r($VOLUMES);
    echo '</pre>';
    
    $fh = fopen("output/volumes.txt", 'w') or die("can't open file");
    foreach($VOLUMES as $value) {

        fwrite($fh, array_search($value, $VOLUMES)."\t".$value."\n");
    }
    fclose($fh);
}
