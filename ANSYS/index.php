<?php 

require_once 'multiGEM.php';

if(!isset($_GET['gemlayers'])) $_GET['gemlayers'] = 1;

if(isset($_POST['submit'])) {

	$gemLayers			= $_POST['gemlayers'];
	$GEM['pitch']		= $_POST['pitch'];
	$GEM['kapton']		= $_POST['kapton'];
	$GEM['metal']		= $_POST['metal'];
	$GEM['outdia']		= $_POST['outdia'];
	$GEM['middia']		= $_POST['middia'];
	$GEM['rim']			= $_POST['rim'];
	$GEM['DRIFTV']		= $_POST['DRIFTV'];
	$GEM['DRIFTT']		= $_POST['DRIFTT'];
	$GEM['INDUCTV']		= $_POST['INDUCTV'];
	$GEM['INDUCTT']		= $_POST['INDUCTT'];

	$totalThickness = $GEM['DRIFTT'] + $GEM['INDUCTT'] + $gemLayers*($GEM['kapton'] + 2*$GEM['metal']);
	$totalVoltage = $GEM['DRIFTV'] + $GEM['INDUCTV'];
	for($i=1; $i<$gemLayers; $i++) {

		$GEM['VGEM'.$i]      = $_POST['VGEM'.$i];
 		$GEM['TRANSFERT'.$i] = $_POST['TRANSFERT'.$i];
		$GEM['TRANSFERV'.$i] = $_POST['TRANSFERV'.$i];
		$totalThickness = $totalThickness + $GEM['TRANSFERT'.$i];
		$totalVoltage = $totalVoltage + $GEM['VGEM'.$i] + $GEM['TRANSFERV'.$i];
	}
	$GEM['VGEM'.$gemLayers] = $_POST['VGEM'.$gemLayers];
	$totalVoltage = $totalVoltage + $GEM['VGEM'.$gemLayers];
	
	// Write configuration to file
	$fh = fopen("output/geometry.txt", 'w') or die("can't open file");
	
    $stringData = "layers ".($_POST['gemlayers'])."\n";
	fwrite($fh, $stringData);
	$stringData = "kapton ".($_POST['kapton']/10)."\n";
	fwrite($fh, $stringData);
	$stringData = "metal ".($_POST['metal']/10)."\n";
	fwrite($fh, $stringData);
	$stringData = "pitch ".($_POST['pitch']/10)."\n";
	fwrite($fh, $stringData);
	$stringData = "outdia ".($_POST['outdia']/10)."\n";
	fwrite($fh, $stringData);
	$stringData = "middia ".($_POST['middia']/10)."\n";
	fwrite($fh, $stringData); 
	$stringData = "rim ".($_POST['rim']/10)."\n";
    fwrite($fh, $stringData);
	$stringData = "driftRegionVoltage ".($_POST['DRIFTV'])."\n";
	fwrite($fh, $stringData);
	$stringData = "driftRegionThickness ".($_POST['DRIFTT']/10)."\n";
	fwrite($fh, $stringData);
	for($i=1; $i<$gemLayers; $i++) {
		$stringData = "GEMVoltage".$i." ".($_POST['VGEM'.$i])."\n";
		fwrite($fh, $stringData);
		$stringData = "transferVoltage".$i." ".($_POST['TRANSFERV'.$i])."\n";
		fwrite($fh, $stringData);
		$stringData = "transferThickness".$i." ".($_POST['TRANSFERT'.$i]/10)."\n";
		fwrite($fh, $stringData);
	}
		$stringData = "GEMVoltage".$gemLayers." ".($_POST['VGEM'.$gemLayers])."\n";
		fwrite($fh, $stringData);
		
	$stringData = "inductRegionVoltage ".($_POST['INDUCTV'])."\n";
	fwrite($fh, $stringData);
	$stringData = "inductRegionThickness ".($_POST['INDUCTT']/10)."\n";
	fwrite($fh, $stringData);
	$stringData = "totalThickness ".($totalThickness/10)."\n";
	fwrite($fh, $stringData);
	$stringData = "totalVoltage ".($totalVoltage)."\n";
	fwrite($fh, $stringData);

	fclose($fh);
}
else {

	$GEM['pitch']	= 0.140;
	$GEM['kapton']	= 0.05;
	$GEM['metal']	= 0.005;
	$GEM['outdia']	= 0.07;
	$GEM['middia']	= 0.05;
	$GEM['rim']		= 0.08;
}

?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1" />
<style type="text/css">
body, html { 
    width: 100%;
    margin: 0px; 
    padding: 0px; 
    border: 0px;
    font-family: arial;
}

input {
    
    border: 1px solid #787878;
}

h1 { padding: 10px; }

</style>
<title>MultiGEM ANSYS configuration file</title>

</head>

<body>

    <h1>MultiGEM ANSYS configuration file generation</h1>
    <hr />
      
    <div style="padding-left: 20px; width: 100%;">
        
        <h3>MultiGEM configuration:</h3>    
        <form action="" method="post">

            <table cellpadding="" cellspacing="">

                <tr>
                    <td width="240px">Gem layers</td>
                    <td>
                        <select name="gemlayers" style="width: 175px" onchange="if (this.value) window.location.href='index.php?gemlayers=' + this.value">
                            <?php 
                            for($i=1; $i<11; $i++) {
                                
                                $sel = ($_GET['gemlayers'] == $i) ? 'selected="selected"' : '';
                                echo '<option '.$sel.' value="'.$i.'">'.$i.'</option>';
                            }
                            ?>  
                        </select>
                    </td>
                </tr>
                <tr>
                    <td>Kapton thickness (K)</td>
                    <td><input type="text" name="kapton" value="<?php echo $GEM['kapton']; ?>" /> mm</td>
                </tr>
                <tr>
                    <td>Metal thickness</td>
                    <td><input type="text" name="metal" value="<?php echo $GEM['metal']; ?>" /> mm</td>
                </tr>
                <tr>
                    <td>Pitch (P)</td>
                    <td><input type="text" name="pitch" value="<?php echo $GEM['pitch']; ?>" /> mm</td>
                </tr>
                <tr>
                    <td>Outer diameter cone (D)</td>
                    <td><input type="text" name="outdia" value="<?php echo $GEM['outdia']; ?>" /> mm</td>
                </tr>
                <tr>
                    <td>Inner diameter cone (d)</td>
                    <td><input type="text" name="middia" value="<?php echo $GEM['middia']; ?>" /> mm</td>
                </tr>
                <tr>
                    <td>Rim</td>
                    <td><input type="text" name="rim" value="<?php echo $GEM['rim']; ?>" /> mm</td>
                </tr>
                <tr>
                    <td>Drift region voltage</td>
                    <td><input type="text" name="DRIFTV" value="<?php echo $GEM['DRIFTV']; ?>" /> V</td>
                </tr>  
                <tr>   
                    <td>Drift region thickness</td>
                    <td><input type="text" name="DRIFTT" value="<?php echo $GEM['DRIFTT']; ?>" /> mm</td>
                </tr>
                <?php
                for($i=1; $i<$_GET['gemlayers']; $i++) {
                    
                    ?>
                    <tr>
                        <td>VGEM <?=$i?></td>
                        <td><input type="text" name="VGEM<?=$i?>" value="<?php echo $GEM['VGEM'.$i]; ?>" /> V</td>
                    </tr>
                    <tr>
                        <td>Transfer voltage <?=$i?></td>
                        <td><input type="text" name="TRANSFERV<?=$i?>" value="<?php echo $GEM['TRANSFERV'.$i]; ?>" /> V</td>
                    </tr>
                    <tr>
                        <td>Transfer thickness <?=$i?></td>
                        <td><input type="text" name="TRANSFERT<?=$i?>" value="<?php echo $GEM['TRANSFERT'.$i]; ?>" /> mm</td>
                    </tr>
                    <?php
                }
                ?>
                <tr>
                    <td>VGEM <?=$_GET['gemlayers']?></td>
                    <td><input type="text" name="VGEM<?=$_GET['gemlayers']?>" value="<?php echo $GEM['VGEM'.$_GET['gemlayers']]; ?>" /> V</td>
                </tr>
                <tr>
                    <td>Induction region voltage</td>
                    <td><input type="text" name="INDUCTV" value="<?php echo $GEM['INDUCTV']; ?>" /> V</td>
                </tr>
                <tr>
                    <td>Induction region thickness</td>
                    <td><input type="text" name="INDUCTT" value="<?php echo $GEM['INDUCTT']; ?>" /> mm</td>
                </tr>
                <tr>
                    <td>&nbsp;</td>
                    <td><input type="submit" name="submit" value="Generate file" /></td>
                </tr>

            </table>

        </form>
        
    </div>

    <br />
    <hr />
    
    <div style="padding-left: 20px; width: 100%;">
    
    <h3>Output code</h3>
    <?php 
    if(isset($_POST['submit'])) {
    

	echo generateFile($GEM, $_POST['gemlayers']);
      	
    }
    ?>
    
    </div>


</body>

</html>
