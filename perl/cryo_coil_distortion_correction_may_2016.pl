#!/usr/bin/perl

use strict;
use Getopt::Tabular;
use MNI::Startup;
use MNI::PathUtilities qw(split_path);
use MNI::Spawn;
MNI::Spawn::SetOptions (verbose => 0);


my $version = 0.1;
my $usage = "\nDistortion correction for images acquired using the cryo coils
on the Bruker scanner (created in May 2016)

$0 [options] -output-dir /path/to/output input_1.mnc input_2.mnc ... input_n.mnc

By default the distortion corrections are sent to the farms (qbatch). If you want to run the corrections on your machine, run:

$0 [options] -spawn -output-dir /path/to/output input_1.mnc input_2.mnc ... input_n.mnc
";

my $help = <<HELP;

Apply distortion correction to cryo-coil files (4 coil array, Bruker scanner)

HELP

Getopt::Tabular::SetHelp($help, $usage);

my $output_dir = undef;
my $resampling = "-sinc";
my $remove_tmp_like = 1;
my $qbatch = 1;
my $spawn = 0;
#my $check_mm_variables = 1;

# handle arguments
my @left_over_args;
my @arg_table = 
    (@DefaultArgs,
      ["-spawn", "boolean", 1, \$spawn,
       "Use the perl system interface to spawn job (run jobs on your local machine) [default]"],
      ["-qbatch", "boolean", 1, \$qbatch,
       "Use qbatch to spawn jobs [default]"],
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
# from the brkr:coil attribute in the MINC header. 

RegisterPrograms(["mincinfo", "mincresample", "qbatch"]);
system("mkdir -p $output_dir") unless (-d $output_dir);

foreach my $mouse (@mice) 
{
	#check whether a coil number for this mouse was specified
	my @parts = split(":", $mouse);
	my $set_coil = undef;
	if($parts[1] == 1 or $parts[1] == 2 or $parts[1] == 3 or $parts[1] == 4)
	{
		$set_coil = $parts[1];
		$mouse = $parts[0];
	}
	#if a coil was specified explicitly, but wasn't on the the 4 coils exit:
	if($#parts == 1 and !$set_coil)
	{
		print "\n\nError: the specified coil number should be one from 1-4. Was $parts[1]. File: ${mouse}\n";
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
	$output .= ".may2016_dist_corr.mnc";
	$output = "${output_dir}/${output}";
	
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
			Spawn(["mincinfo", "-attval", "brkr:coil", $mouse],  stdout=>\${coil});
		};
		if($@)
		{
			print "\nError: the brkr:coil entry for file $mouse must be one from 1-4 (was: $coil )\n\n";
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
      ${coil} == 4 ) {
    $trans_to_use = "/hpf/largeprojects/MICe/matthijs/2016-05-cryo-distortion-correction/correction_coil${coil}/coil${coil}_lsq9_scaling_only.xfm";
  }
  else {
    print "\nError: the brkr:coil entry for file $mouse must be one from 1-4 (was: $coil )\n";
    next;
  }
  $correction_command .= " -transform ${trans_to_use} ";
  
  
  
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
    $combined_call .= " \" | qbatch - --mem=2G -N dc-cryo-2G-${base}";
  }
	

	
	if($Execute) 
	{
      system("$combined_call");
 	}
  else 
  {
    print "$combined_call\n";
  }
}

print("\n");

exit;

