#!/usr/bin/perl

use strict;
use Getopt::Tabular;
use MNI::Startup;
use MNI::PathUtilities qw(split_path);
use MNI::Spawn;
MNI::Spawn::SetOptions (verbose => 0);


my $version = 0.1;
my $usage = "\nDistortion correction for images acquired using the gradient we
received at MICe in the fall of 2013. 
$0 [options] -output-dir /path/to/output input_1.mnc input_2.mnc ... input_n.mnc

By default the distortion corrections are sent to the farms (qbatch). If you want to run the corrections on your machine, run:

$0 [options] -spawn -output-dir /path/to/output input_1.mnc input_2.mnc ... input_n.mnc
";

my $help = <<HELP;

Apply distortion correction to brain-in-the-basket files (16 coil array, fall 2013 and onwards)

HELP

Getopt::Tabular::SetHelp($help, $usage);

my $output_dir = undef;
my $resampling = "-sinc";
my $use_lambda_0_01_distance_10 = 0;
my $remove_tmp_like = 1;
my $qbatch = 1;
my $spawn = 0;

# handle arguments
my @left_over_args;
my @arg_table = 
    (@DefaultArgs,
      ["-spawn", "boolean", 1, \$spawn,
       "Use the perl system interface to spawn job (run jobs on your local machine)"],
      ["-qbatch", "boolean", 1, \$qbatch,
       "Use qbatch to spawn jobs [default]"],
      ["-sinc", "const", "-sinc", \$resampling,
       "Use sinc resampling [default]"],
      ["-trilinear", "const", "-trilinear", \$resampling,
       "Use trilinear resampling"],
      ["-tricubic", "const", "-tricubic", \$resampling,
       "Use tricubic interpolation"],
      ["-output-dir", "string", 1, \$output_dir,
       "Directory in which to place the output files."],
#       ["-remove-temp-likefile|-no-remove-temp-likefile", "boolean", 1, \$remove_tmp_like,
#        "Remove the temporary like file that is created for the resampling."],
      );
GetOptions(\@arg_table, \@ARGV, \@left_over_args) or die "\n";

my @mice = @left_over_args;

#if the -output-dir was not used, but the first argument
#is a directory, use this value to get the output directory
if(!$output_dir and -d $mice[0])
{
	$output_dir = shift @mice;
}

if(!$output_dir)
{
	print "\n\nError: please specify an output directory.\n\n";
	die $usage;
}

die $usage unless ($output_dir and $#mice >= 0);



# 1) Perform the correct transformation to the files
#
# The coil each file was scanned in is retrieved
# from the vnmr:coil attribute in the MINC header. 

RegisterPrograms(["mincinfo", "mincresample", "qbatch"]);
system("mkdir -p $output_dir") unless (-d $output_dir);

foreach my $mouse (@mice) 
{
	#check whether a coil number for this mouse was specified
	my @parts = split(":", $mouse);
	my $set_coil = undef;
	if($parts[1] == 1 or $parts[1] == 2 or $parts[1] == 3 or $parts[1] == 4 or $parts[1] == 5 or
	   $parts[1] == 6 or $parts[1] == 7 or $parts[1] == 8 or $parts[1] == 9 or $parts[1] == 10 or
	   $parts[1] == 11 or $parts[1] == 12 or $parts[1] == 13 or $parts[1] == 14 or 
	   $parts[1] == 15 or $parts[1] == 16)
	{
		$set_coil = $parts[1];
		$mouse = $parts[0];
	}
	#if a coil was specified explicitly, but wasn't on the the 14 working coils exit:
	if($#parts == 1 and !$set_coil)
	{
		print "\n\nError: the specified coil number should be one from 1-16. Was $parts[1]. File: ${mouse}\n";
		next;
	}
		
	my ($dir, $base, $ext) = split_path($mouse, 'last', [qw(gz z Z)]);	
	#at this point, $mouse should hold a minc file, verify this:
	if(!($ext eq '.mnc'))
	{
		print "\n\nError: $mouse is not a minc file\n";
		next;
	}
	my $output = $base;
	$output .= ".sept2014_dist_corr.mnc";
	$output = "${output_dir}/${output}";
	
  # the deformation grid, needed to find out whether the input file actually falls
  # within the correction area
  my $deformation_grid;
	
	#transform options
	my $correction_command = "mincresample ${resampling} -2 ";
  if($Clobber) {
    $correction_command .= " -clobber ";
  }
	
	my $coil;
	
	my $break_to_next = 0;
	if(!$set_coil)
	{
		eval
		{
			Spawn(["mincinfo", "-attval", "vnmr:coil", $mouse],  stdout=>\${coil});
		};
		if($@)
		{
			print "\nError: the vnmr:coil entry for file $mouse must be one from 1-16 (was: $coil )\n\n";
			$break_to_next = 1;
		}
		$coil =~ s/\n//;
		$coil =~ s/ //g;
		print "\nFound coil:   $coil    for file: $mouse    ";
	}
	else
	{
		$coil = $set_coil
	}
	
	# go to next in loop 
	if($break_to_next == 1){
    next;
  }

  ##############################################################################
  # Deal with the transformation
  my $trans_to_use;
  #
  # For coils 1-4 and 13-16, we use different settings for tagtoxfm_bspline
  # because the tagpoints chosen for these coils lie further apart, less
  # smoothing is necessary when approximating the bsplines
  #
  # Final settings are: lambda 0.001, distance 10, range [-5,5], tags manually placed by
  #                     Jeff Emmanuel, and a large likefile
  #                     for the deformation field
  #
  #
  if( ${coil} == 1 or
      ${coil} == 2 or 
      ${coil} == 3 or 
      ${coil} == 4 or 
      ${coil} == 13 or 
      ${coil} == 14 or 
      ${coil} == 15 or 
      ${coil} == 16 ) {
    $trans_to_use = "/hpf/largeprojects/MICe/tools/distortion_correction/2013-11-Distortion-Correction/2014-09-bsplines-on-landmarks/CONCAT_CT_and_MR_coil_${coil}_lambda_0.001_distance_10_range_5.xfm";
    $deformation_grid = "/hpf/largeprojects/MICe/tools/distortion_correction/2013-11-Distortion-Correction/2014-09-bsplines-on-landmarks/CONCAT_CT_and_MR_coil_${coil}_lambda_0.001_distance_10_range_5_grid_0.mnc";
  }
  #
  # Coils 8-12
  #
  # Final settings are: lambda 0.01, distance 10, range [-5,5], tags manually place by
  #                     Jeff Emmanuel, and a large likefile
  #                     for the deformation field
  #
  #
  elsif (${coil} == 5 or
        ${coil} == 6 or 
        ${coil} == 7 or 
        ${coil} == 8 or 
        ${coil} == 9 or 
        ${coil} == 10 or 
        ${coil} == 11 or 
        ${coil} == 12 ) { 
    $trans_to_use = "/hpf/largeprojects/MICe/tools/distortion_correction/2013-11-Distortion-Correction/2014-09-bsplines-on-landmarks/CONCAT_CT_and_MR_coil_${coil}_lambda_0.01_distance_10_range_5.xfm";
    $deformation_grid = "/hpf/largeprojects/MICe/tools/distortion_correction/2013-11-Distortion-Correction/2014-09-bsplines-on-landmarks/CONCAT_CT_and_MR_coil_${coil}_lambda_0.01_distance_10_range_5_grid_0.mnc";
  }
  else {
    print "\nError: the vnmr:coil entry for file $mouse must be one from 1-16 (was: $coil )\n";
    next;
  }
  $correction_command .= " -transform ${trans_to_use} ";
  
  # check the overlap between the bounding box of the disortion 
  #
  # Transformation
  #
  my ($trans_x_start, $trans_y_start, $trans_z_start);
  Spawn(["mincinfo", "-attval", "xspace:start", $deformation_grid],  stdout=>\${trans_x_start});
  Spawn(["mincinfo", "-attval", "yspace:start", $deformation_grid],  stdout=>\${trans_y_start});
  Spawn(["mincinfo", "-attval", "zspace:start", $deformation_grid],  stdout=>\${trans_z_start});
  my ($trans_x_step, $trans_y_step, $trans_z_step);
  Spawn(["mincinfo", "-attval", "xspace:step", $deformation_grid],  stdout=>\${trans_x_step});
  Spawn(["mincinfo", "-attval", "yspace:step", $deformation_grid],  stdout=>\${trans_y_step});
  Spawn(["mincinfo", "-attval", "zspace:step", $deformation_grid],  stdout=>\${trans_z_step});
  my ($trans_x_length, $trans_y_length, $trans_z_length);
  Spawn(["mincinfo", "-attval", "xspace:length", $deformation_grid],  stdout=>\${trans_x_length});
  Spawn(["mincinfo", "-attval", "yspace:length", $deformation_grid],  stdout=>\${trans_y_length});
  Spawn(["mincinfo", "-attval", "zspace:length", $deformation_grid],  stdout=>\${trans_z_length});
  my ($trans_x_min, $trans_x_max, $trans_y_min, $trans_y_max, $trans_z_min, $trans_z_max);
  if($trans_x_step < 0) {
    $trans_x_min = $trans_x_start + ($trans_x_step * $trans_x_length);
    $trans_x_max = $trans_x_start + 0;
  }
  else {
    $trans_x_min = $trans_x_start + 0;
    $trans_x_max = $trans_x_start + ($trans_x_step * $trans_x_length);
  }
  if($trans_y_step < 0) {
    $trans_y_min = $trans_y_start + ($trans_y_step * $trans_y_length);
    $trans_y_max = $trans_y_start + 0;
  }
  else {
    $trans_y_min = $trans_y_start + 0;
    $trans_y_max = $trans_y_start + ($trans_y_step * $trans_y_length);
  }
  if($trans_z_step < 0) {
    $trans_z_min = $trans_z_start + ($trans_z_step * $trans_z_length);
    $trans_z_max = $trans_z_start + 0;
  }
  else {
    $trans_z_min = $trans_z_start + 0;
    $trans_z_max = $trans_z_start + ($trans_z_step * $trans_z_length);
  }
  
  #
  # Input file
  #
  my ($input_x_start, $input_y_start, $input_z_start);
  Spawn(["mincinfo", "-attval", "xspace:start", $mouse],  stdout=>\${input_x_start});
  Spawn(["mincinfo", "-attval", "yspace:start", $mouse],  stdout=>\${input_y_start});
  Spawn(["mincinfo", "-attval", "zspace:start", $mouse],  stdout=>\${input_z_start});
  my ($input_x_step, $input_y_step, $input_z_step);
  Spawn(["mincinfo", "-attval", "xspace:step", $mouse],  stdout=>\${input_x_step});
  Spawn(["mincinfo", "-attval", "yspace:step", $mouse],  stdout=>\${input_y_step});
  Spawn(["mincinfo", "-attval", "zspace:step", $mouse],  stdout=>\${input_z_step});
  my ($input_x_length, $input_y_length, $input_z_length);
  Spawn(["mincinfo", "-attval", "xspace:length", $mouse],  stdout=>\${input_x_length});
  Spawn(["mincinfo", "-attval", "yspace:length", $mouse],  stdout=>\${input_y_length});
  Spawn(["mincinfo", "-attval", "zspace:length", $mouse],  stdout=>\${input_z_length});
  my ($input_x_min, $input_x_max, $input_y_min, $input_y_max, $input_z_min, $input_z_max);
  if($input_x_step < 0) {
    $input_x_min = $input_x_start + ($input_x_step * $input_x_length);
    $input_x_max = $input_x_start + 0;
  }
  else {
    $input_x_min = $input_x_start + 0;
    $input_x_max = $input_x_start + ($input_x_step * $input_x_length);
  }
  if($input_y_step < 0) {
    $input_y_min = $input_y_start + ($input_y_step * $input_y_length);
    $input_y_max = $input_y_start + 0;
  }
  else {
    $input_y_min = $input_y_start + 0;
    $input_y_max = $input_y_start + ($input_y_step * $input_y_length);
  }
  if($input_z_step < 0) {
    $input_z_min = $input_z_start + ($input_z_step * $input_z_length);
    $input_z_max = $input_z_start + 0;
  }
  else {
    $input_z_min = $input_z_start + 0;
    $input_z_max = $input_z_start + ($input_z_step * $input_z_length);
  }
  
  #
  # Most basic test: check extreme limits only
  #
  if($input_x_min > $trans_x_max or $input_x_max < $trans_x_min) {
    print("\nError: the x-dimension of the input file does not fall into the correction field: ", $mouse, "\n");
    next;
  }
  if($input_y_min > $trans_y_max or $input_y_max < $trans_y_min) {
    print("\nError: the y-dimension of the input file does not fall into the correction field: ", $mouse, "\n");
    next;
  }
  if($input_z_min > $trans_z_max or $input_z_max < $trans_z_min) {
    print("\nError: the z-dimension of the input file does not fall into the correction field: ", $mouse, "\n");
    next;
  }
  
  
  # the corrections were made in native coil space, so we don't need any other like file than the original file:
  $correction_command .= " -like ${mouse} ";
  $correction_command .= " ${mouse} ${output} ";

  my $full_line_of_commands = "${correction_command}";
  
  
  # add information about the distortion correction:
  my $minc_header_addon = "";
  if($qbatch and !$spawn)  {
    $minc_header_addon .= "minc_modify_header -sinsert distortioncorrection:coil=${coil} -sinsert distortioncorrection:programcall=\\\"";
  }
  else{
    $minc_header_addon .= "minc_modify_header -sinsert distortioncorrection:coil=${coil} -sinsert distortioncorrection:programcall=\"";
  }
  $minc_header_addon .= $0;
  $minc_header_addon .= " ";
  for(my $i = 0; $i <= $#ARGV; $i++){
    $minc_header_addon .= $ARGV[$i];
    $minc_header_addon .= " ";
  }
  if($qbatch and !$spawn)  {
    $minc_header_addon .= "\\\" ${output} ";
  }
  else{
    $minc_header_addon .= "\" ${output} ";
  }
	

	my $combined_call = "";
	if($qbatch and !$spawn)  {
    $combined_call .= " echo \""
}
	$combined_call .= $full_line_of_commands;
	$combined_call .= ";";
	$combined_call .= $minc_header_addon;
	if($qbatch and !$spawn)  {
    $combined_call .= " \" | qbatch - --mem=4G -N dc-4G-${base}";
  }


	
	if($Execute) 
	{
      system("$combined_call")
 	}
  else 
  {
    print "$combined_call\n";
  }
}

print("\n");

exit;

