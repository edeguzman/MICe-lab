#!/usr/bin/perl

use strict;
use Getopt::Tabular;
use MNI::Startup;
use MNI::PathUtilities qw(split_path);
use MNI::Spawn;
MNI::Spawn::SetOptions (verbose => 0);


my $version = 0.1;
my $usage = "\nDistortion correction for images acquired using the saddle coils
on the Varian scanner (created in August 2015)

$0 [options] -output-dir /path/to/output input_1.mnc input_2.mnc ... input_n.mnc

By default the distortion corrections are sent to the farms (sge). If you want to run the corrections on your machine, run:

$0 [options] -spawn -output-dir /path/to/output input_1.mnc input_2.mnc ... input_n.mnc
";

my $help = <<HELP;

Apply distortion correction to saddle-coil files (7 coil array, Varian scanner)

HELP

Getopt::Tabular::SetHelp($help, $usage);

my $output_dir = undef;
my $resampling = "-sinc";
my $remove_tmp_like = 1;
my $sge = 1;
my $spawn = 0;
my $sgequeue = "all.q";
my $check_mm_variables = 1;

# handle arguments
my @left_over_args;
my @arg_table = 
    (@DefaultArgs,
      ["-spawn", "boolean", 1, \$spawn,
       "Use the perl system interface to spawn job (run jobs on your local machine)"],
      ["-sge", "boolean", 1, \$sge,
       "Use SGE to spawn jobs [default]"],
      ["-check_mm|-no-check_mm", "boolean", 1, \$check_mm_variables,
       "Use the mm_ppe2, mm_ppe and mm_pro header values to determine that the correct hive table was used in the file reconstruction"],
      ["-sge-queue", "string", 1, \$sgequeue,
       "Which SGE queues to send the jobs to [default: all.q]"],
      ["-sinc", "const", "-sinc", \$resampling,
       "Use sinc resampling [default]"],
      ["-trilinear", "const", "-trilinear", \$resampling,
       "Use trilinear resampling"],
      ["-tricubic", "const", "-tricubic", \$resampling,
       "Use tricubic interpolation"],
      ["-output-dir", "string", 1, \$output_dir,
       "Directory in which to place the output files."]
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

RegisterPrograms(["mincinfo", "mincresample", "sge_batch"]);
system("mkdir -p $output_dir") unless (-d $output_dir);

foreach my $mouse (@mice) 
{
	#check whether a coil number for this mouse was specified
	my @parts = split(":", $mouse);
	my $set_coil = undef;
	if($parts[1] == 1 or $parts[1] == 2 or $parts[1] == 3 or $parts[1] == 4 or $parts[1] == 5 or
	   $parts[1] == 6 or $parts[1] == 7)
	{
		$set_coil = $parts[1];
		$mouse = $parts[0];
	}
	#if a coil was specified explicitly, but wasn't on the the 14 working coils exit:
	if($#parts == 1 and !$set_coil)
	{
		print "\n\nError: the specified coil number should be one from 1-7. Was $parts[1]. File: ${mouse}\n";
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
	$output .= ".aug2015_dist_corr.mnc";
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
			print "\nError: the vnmr:coil entry for file $mouse must be one from 1-7 (was: $coil )\n\n";
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
  if (${coil} == 1 or
      ${coil} == 2 or 
      ${coil} == 3 or 
      ${coil} == 4 or
      ${coil} == 5 or
      ${coil} == 6 or 
      ${coil} == 7 ) {
    $trans_to_use     = "/projects/muisjes/backup/matthijs/2015-07-saddle-distortion-correction/final_corrections/distortion_correction_coil_${coil}_lambda_0.01_distance_10.xfm";
    $deformation_grid = "/projects/muisjes/backup/matthijs/2015-07-saddle-distortion-correction/final_corrections/distortion_correction_coil_${coil}_lambda_0.01_distance_10_grid_0.mnc";
  }
  else {
    print "\nError: the vnmr:coil entry for file $mouse must be one from 1-7 (was: $coil )\n";
    next;
  }
  $correction_command .= " -transform ${trans_to_use} ";
  
  
  # use header information about the file's reconstruction to check that the correct hive table was used
  if ($check_mm_variables) {
    my ($file_mm_ppe2, $file_mm_ppe, $file_mm_pro);
    Spawn(["mincinfo", "-attval", "vnmr:mm_ppe2", $mouse], stdout=>\${file_mm_ppe2});
    Spawn(["mincinfo", "-attval", "vnmr:mm_ppe", $mouse],  stdout=>\${file_mm_ppe});
    Spawn(["mincinfo", "-attval", "vnmr:mm_pro", $mouse],  stdout=>\${file_mm_pro});
    if ($coil == 1) {
      if($file_mm_ppe2 != -0.0 or
         $file_mm_ppe  !=  0.0 or 
         $file_mm_pro  !=  2.5){
        print "\nError: mm_ values for file $mouse incorrect. For coil $coil this should be: (mm_ppe2 = -0.0, mm_ppe = 0.0, mm_pro = 2.5), was: ${file_mm_ppe2}, ${file_mm_ppe}, ${file_mm_pro}\n";
        next;
      }
    }
    if ($coil == 2) {
      if($file_mm_ppe2 != -2.6 or
         $file_mm_ppe  !=  5.0 or 
         $file_mm_pro  != -2.5){
        print "\nError: mm_ values for file $mouse incorrect. For coil $coil this should be: (mm_ppe2 = -2.6, mm_ppe = 5.0, mm_pro = -2.5), was: ${file_mm_ppe2}, ${file_mm_ppe}, ${file_mm_pro}\n";
        next;
      }
    }
    if ($coil == 3) {
      if($file_mm_ppe2 != -5.4 or
         $file_mm_ppe  != -0.0 or 
         $file_mm_pro  !=  2.5){
        print "\nError: mm_ values for file $mouse incorrect. For coil $coil this should be: (mm_ppe2 = -5.4, mm_ppe = -0.0, mm_pro = 2.5), was: ${file_mm_ppe2}, ${file_mm_ppe}, ${file_mm_pro}\n";
        next;
      }
    }
    if ($coil == 4) {
      if($file_mm_ppe2 != -2.9 or
         $file_mm_ppe  != -5.0 or 
         $file_mm_pro  != -2.5){
        print "\nError: mm_ values for file $mouse incorrect. For coil $coil this should be: (mm_ppe2 = -2.9, mm_ppe = -5.0, mm_pro = -2.5), was: ${file_mm_ppe2}, ${file_mm_ppe}, ${file_mm_pro}\n";
        next;
      }
    }
    if ($coil == 5) {
      if($file_mm_ppe2 !=  2.7 or
         $file_mm_ppe  != -5.0 or 
         $file_mm_pro  != -2.5){
        print "\nError: mm_ values for file $mouse incorrect. For coil $coil this should be: (mm_ppe2 = 2.7, mm_ppe = -5.0, mm_pro = -2.5), was: ${file_mm_ppe2}, ${file_mm_ppe}, ${file_mm_pro}\n";
        next;
      }
    }
    if ($coil == 6) {
      if($file_mm_ppe2 !=  5.55 or
         $file_mm_ppe  != -0.1 or 
         $file_mm_pro  !=  2.5){
        print "\nError: mm_ values for file $mouse incorrect. For coil $coil this should be: (mm_ppe2 = 5.55, mm_ppe = -0.1, mm_pro = 2.5), was: ${file_mm_ppe2}, ${file_mm_ppe}, ${file_mm_pro}\n";
        next;
      }
    }
    if ($coil == 7) {
      if($file_mm_ppe2 !=  2.8 or
         $file_mm_ppe  !=  4.8 or 
         $file_mm_pro  != -2.5){
        print "\nError: mm_ values for file $mouse incorrect. For coil $coil this should be: (mm_ppe2 = 2.8, mm_ppe = 4.8, mm_pro = -2.5), was: ${file_mm_ppe2}, ${file_mm_ppe}, ${file_mm_pro}\n";
        next;
      }
    }
  }
  
  
  
  # check the overlap between the bounding box of the distortion 
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
  # the z dimension is ridiculously large in these files, so we're only checking x and y
#   if($input_z_min > $trans_z_max or $input_z_max < $trans_z_min) {
#     print("\nError: the z-dimension of the input file does not fall into the correction field: ", $mouse, "\n");
#     next;
#   }
  
  
  # the corrections were made in native coil space, so we don't need any other like file than the original file:
  $correction_command .= " -like ${mouse} ";
  $correction_command .= " ${mouse} ${output} ";

  my $full_line_of_commands = "${correction_command}";
  
  
  # add information about the distortion correction:
  my $minc_header_addon = "";
  if($sge and !$spawn)  {
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
  if($sge and !$spawn)  {
    $minc_header_addon .= "\\\" ${output} ";
  }
  else{
    $minc_header_addon .= "\" ${output} ";
  }
	
	my $combined_call = "";
	if($sge and !$spawn)  {
    $combined_call .= "sge_batch -q $sgequeue -l vf=2G -J dc-saddle-2G-${base} \"";
  }
	$combined_call .= $full_line_of_commands;
	$combined_call .= ";";
	$combined_call .= $minc_header_addon;
	if($sge and !$spawn)  {
    $combined_call .= " \"";
  }
	

	
	if($Execute) 
	{
    if(!$sge or $spawn) 
    {
      system("$combined_call")
    }
    else 
    {
      #Spawn(["sge_batch", "-q", $sgequeue, "-l", "vf=4G", "-J", "dc-4G-${base}", $combined_call]);
      system("$combined_call");
      #system("sge_batch -l vf=4G -J dc-4G-${base}  ", "$full_line_of_commands");
    }
 	}
  else 
  {
    print "$combined_call\n";
  }
}

print("\n");

exit;

