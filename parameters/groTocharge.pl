#!/usr/local/bin/perl 
#
##  usage perl program.pl input.gro  
##  Assign integer charges  for Lysine (+1), Argenine (+1), Aspartate (-1), Glutamate (-1)  and  rest others are zero
#
#
if ($#ARGV < 0 ){
    print "usage: perl program.pl protein.gro #gro file has the list of residue and atoms \n";
        exit;
        }
$input_file  = $ARGV[0];

open(ReadGroFile,$input_file) || die "Can't Open File: $input_file.\n";
$line_input = <ReadGroFile>;
$line_input = <ReadGroFile>;
#print "$line_input";
$line_input = <ReadGroFile>;
$charge_count = 0;
while ($line_input ne "")
   {
     chop ($line_input);
     @split_line = split(" ", $line_input);
     $residueID = @split_line[0];
     $residue = @split_line[1];
     $atom  = @split_line[2];
     if($residue =~ /ARG/ && $atom =~ /CB/)
        {
        $charge = 1.0;
        #printf " %6d   %3s  %2s  %8.4f \n", $residueID, $residue, $atom, $charge ;
        @residueID_store[$charge_count] = $residueID;
        @charge_store[$charge_count] = $charge;
        $charge_count++;
        }
     elsif($residue =~ /LYS/ && $atom =~ /CB/)
        {
        $charge = 1.0;
        #printf " %6d   %3s  %2s  %8.4f \n", $residueID, $residue, $atom, $charge ;
        @residueID_store[$charge_count] = $residueID;
        @charge_store[$charge_count] = $charge;
        $charge_count++;
        }
     elsif($residue =~ /ASP/ && $atom =~ /CB/)
        {
        $charge = -1.0;
        #printf " %6d   %3s  %2s  %8.4f \n", $residueID, $residue, $atom, $charge ;
        @residueID_store[$charge_count] = $residueID;
        @charge_store[$charge_count] = $charge;
        $charge_count++;
        }
     elsif($residue =~ /GLU/ && $atom =~ /CB/)
        {
        $charge = -1.0;
        #printf " %6d   %3s  %2s  %8.4f \n", $residueID, $residue, $atom, $charge ;
        @residueID_store[$charge_count] = $residueID;
        @charge_store[$charge_count] = $charge;
        $charge_count++;
        }

     $line_input = <ReadGroFile>;
   }

print "$charge_count \n";
for($i = 0; $i < $charge_count; $i++)
{

printf "%6d   %8.4f \n", @residueID_store[$i], @charge_store[$i];

}

