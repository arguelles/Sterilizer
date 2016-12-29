<?php
	if(!isset($_POST) || !array_key_exists("store",$_POST)
	   || !array_key_exists("dataset",$_POST)
	   || !array_key_exists("name",$_POST) || !array_key_exists("data",$_POST)){
		header('HTTP/1.1 400 Bad Request', true, 400);
		return;
	}
	$datasetDir=$_POST["dataset"];
	if(!is_dir($datasetDir)){
		header('HTTP/1.1 400 Bad Request', true, 400);
		echo "Bad dataset name: ".$datasetDir;
		return;
	}
	$dataFile=$_POST["name"];
	$outputFile=$datasetDir."/".$dataFile;
	
	if(is_file($outputFile)){
		header('HTTP/1.1 409 Conflict', true, 409);
		echo "Data entry ".$datasetDir."/".$dataFile." already exists";
		return;
	}
	$data=$_POST["data"];
	$maxDataSize=1024*10;
	if(strlen($data)>$maxDataSize){
		header('HTTP/1.1 413 Request Entity Too Large', true, 413);
		echo "Data is larger than ".$maxDataSize." bytes";
		return;
	}
	
	file_put_contents($outputFile,$data);
	chmod($outputFile, 0600); //private and NOT EXECUTABLE
	header('HTTP/1.1 201 Created', true, 201);
	echo "Data entry ".$datasetDir."/".$dataFile." was created";
?>