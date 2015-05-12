##############################################################################
# 
# upe.tcl - 
#
# Copyright (C) 2001-2004 Institute of Physical Biology,
#                    Heinrich-Heine Universität Düsseldorf, Germany
#                    <construct@biophys.uni-duesseldorf.de>
#
# This file is part of ConStruct.
# 
# ConStruct is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# ConStruct is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with ConStruct; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# 
###############################################################################


#
#  CVS $Id: upe.tcl,v 1.6 2004-08-11 08:14:22 wilm Exp $    
#

###############################################################################
#               
# main function: upe
#               
#   Calculation of the unbiased probability estimator
#                    =        =           =
#   according to Chiu & Kolodziejczak (1991) CABIOS 7, 347-352.
#       Inferring consensus structure from nucleic acid sequences
#
#
#
namespace eval Upe {
  
  
    variable debug_level   0
    variable verbose_level 0
    
    variable chi2
    variable chi2_p

    
    

   ###   upe   
   #
   # Calculation of the unbiased probability estimator
   #                    =        =           =
   #    according to Chiu & Kolodziejczak (1991) CABIOS 7, 347-352.
   #                 Inferring consensus structure from nucleic acid sequences
   #
   #	For a general introduction to information theory see:
   #		Shannon, C.E. & Weaver, W. (1948).
   #		The mathematical theory of communication.
   #		The University of Illinois Press, Urbana.
   #
   # Returns a list with elements: chi2-prob X:Y
   #    with X:Y the significant nt combinations
   #    if chi2-prob == 0. then there is no significant X:Y
   #
   # Input: n = number of sequences
   #        i = numbers of nts at position i             (this is array fi in calling routine)
   #        j = numbers of nts at position j             (this is array fj in calling routine)
   #        k = numbers of nt pairs at positions i and j (this is array f  in calling routine)
   #            n, i, j, and k might be integers or reals
   #        listNt = list with allowed nts
   #        debug = 0 || 1
   #                1 = true => check of consistency of arrays fi, fj, and f;
   #                            output tables with frequency of nts at pos. i and j,
   #                                               probability of pairs at i:j,
   #                                               interdependence entropies of pairs at i:j (==Shannon's entropy H),
   #                                   the expected mutual information or rate of information transmission I,
   #                                   its significance level, and
   #                                   the significant pairs.
   #        bit = 0 || 1
   #              1 = true  => do calculations with log_2 (see Schneider, Stormo, Gold & Ehrenfeucht (1986))
   #              0 = false => do calculations with log_e (see Chiu & Kolodziejczak (1991))
   #        unbiased = 0 || 1
   #              1 = true  => do calculations using unbiased probability estimation  (see Chiu & Kolodziejczak (1991))
   #                           Approbriate for smaller numbers of sequences and/or for positions where some nt types are not observed
   #              0 = false => do calculations using maximum likelyhood estimation    (see Gutell, Power, Hertz, Putz & Stormo (1992))
   #
   # Globals: chi2 = "table" of chi^2 probabilities; setup by proc SetupChi2Distr;
   #                 used for determination of significance levels;
   #                 with listNt={    A U G C} only the row with  9 degrees of freedom is used;
   #                 with listNt={-   A U G C} only the row with 16 degrees of freedom is used;
   #                 with listNt={- N A U G C} only the row with 25 degrees of freedom is used
   #
   proc upe {n i j k listNt bit unbiased} {
   ########
        variable chi2
        variable chi2_p
    	
        variable debug_level
        variable verbose_level
   
    	
       	upvar $i fi
       	upvar $j fj
       	upvar $k f

       
       	#vputs "upe $x $y $n $i $j $k $listNt $bit $unbiased" 4
       	
        if {$verbose_level>1} {;		# check for correctness/plausibility of values
            CheckPlausibility $listNt fi fj $n f
        }
    		
    
       
        # Calculate single nucleotide dependencies
       	set lenListNt [llength $listNt]
       	if {$unbiased} {
       		set n16 [expr {1.*($n+$lenListNt*$lenListNt)}]
       		set n4  [expr {1.*($n+$lenListNt)}]
       		set One 1.
       	} else {
       		set n16 [expr {1.*$n}]
       		set n4  [expr {1.*$n}]
       		set One 0.
       	}
       	set numSignificantPair 0
       	set significantPair($numSignificantPair) [list 0. a a]
       	if {$verbose_level>1} {
       		puts "   Frequency at"
       		puts "     position"
       		set title "     i      j  "
       		foreach a $listNt {
       			append title " p(x,$a)"
       		}
       		puts "$title"
       	}
       	foreach a $listNt {
       		set pi($a) [expr {($fi($a)+$One)/$n4}]
       		set pj($a) [expr {($fj($a)+$One)/$n4}]
       		#vputs [format "%s %6.3f %6.3f " $a $pi($a) $pj($a)] 7 1
       		foreach b $listNt {
       			set p($a,$b) [expr {($f($a,$b)+$One)/$n16}]
       			#vputs [format "%6.3f " $p($a,$b)] 7
       			if { $p($a,$b)>[lindex $significantPair($numSignificantPair) 0] } {
       				set numSignificantPair 0
       				set significantPair($numSignificantPair) [list $p($a,$b) $a $b]
       			} elseif { $p($a,$b)==[lindex $significantPair($numSignificantPair) 0] } {
       				incr numSignificantPair
       				set significantPair($numSignificantPair) [list $p($a,$b) $a $b]
       			}
       		}
       		#vputs " " 7
       	}
       
       	if {$bit} {;	# do calculations with log_2
       		set loge2 [expr {log(2)}]
       	} else {;	# do calculations with log_e
       		set loge2 1.
       	}
       
       	set Hi 0.
       	set Hj 0.
       	foreach a $listNt {
       		set pi($a) [expr {($fi($a)+$One)/$n4}]
       		if {$pi($a)>0.} {
       			set logpi($a) [expr {log($pi($a))}]
       			set Hi [expr {$Hi - $pi($a)*$logpi($a)/$loge2}]
       		} else {
       			set logpi($a) 0.
       		}
       		set pj($a) [expr {($fj($a)+$One)/$n4}]
       		if {$pj($a)>0.} {
       			set logpj($a) [expr {log($pj($a))}]
       			set Hj [expr {$Hj - $pj($a)*$logpj($a)/$loge2}]
       		} else {
       			set logpj($a) 0.
       		}
       	}
       	#vputs "  =============" 7
       	#vputs [format "   %.3f  %.3f  H(i/j)" $Hi $Hj] 7
       	#vputs " " 7
       
       	set sum 0.
       	if {$verbose_level>1} {
       		set title " "
       		foreach a $listNt {
       			append title " H(x,$a)"
       		}
       		append title "  sum"
       		puts "$title"
       	}
       	
       	foreach a $listNt {
       		#vputs [format "%s " $a] 7 1
       		set sum1 0.
       		foreach b $listNt {;	# for "biased" estimation p=0 is allowed
       			if {$p($a,$b)>0. && $pi($a)>0. && $pj($b)>0.} {
       			set sum0 [expr {$p($a,$b)*(log($p($a,$b))-$logpi($a)-$logpj($b))/$loge2}]
       			# set sum0 [expr $p($a,$b)*log($p($a,$b)/($pi($a)*$pj($b)))/$loge2]; # this is equivalent to the above line, but slower
       			} else {
       				set sum0 0.
       			}
       			#vputs [format "%6.3f " $sum0] 7 1
       			set sum1  [expr {$sum1 + $sum0}]
       		}
       		#vputs [format "%6.3f" $sum1] 7
       		set sum [expr {$sum + $sum1}]
       	}
       	#vputs [format "%[expr {$lenListNt*7+2+6}]s" ======] 7
       
       	#vputs [format "%[expr {$lenListNt*7-3}]s I = %6.3f" " " $sum ] 7 1
       	lappend retlist $sum
       	set df [expr {$lenListNt*$lenListNt - 2*$lenListNt +1}]
       	set signifikant 1
       	foreach prob $chi2_p {
       		if {$sum>=[expr {$chi2($df,$prob)/(2.*$n)}]} {
       			#vputs [format " >= %.3f = x2(%2d,$prob) !!" [expr {$chi2($df,$prob)/(2.*$n)}] $df] 7
       			lappend retlist [format " >= %.3f = x2(%2d,$prob)"    [expr {$chi2($df,$prob)/(2.*$n)}] $df]
       			set signifikant 0
       			break
       		}
       	}
       	if {$signifikant==1} {
       		#vputs " => not significant" 7
       		lappend retlist " => not significant"
       	}
       
       	if {$verbose_level>1} {
       		if {$Hi>0.} {
       			puts [format "%[expr {$lenListNt*7-8}]s R1(ij) = %6.3f" " " [expr {$sum/$Hi}] ]
       		} else {
       			puts [format "%[expr {$lenListNt*7-8}]s R1(ij) = ???"   " "                 ]
       		}
       		if {$Hj>0.} {
       			puts [format "%[expr {$lenListNt*7-8}]s R2(ij) = %6.3f" " " [expr {$sum/$Hj}] ]
       		} else {
       			puts [format "%[expr {$lenListNt*7-8}]s R2(ij) = ???"   " "                 ]
       		}
       	}
       	if {$Hi>0.} {lappend retlist [expr {$sum/$Hi}]} else {lappend retlist 0.}
       	if {$Hj>0.} {lappend retlist [expr {$sum/$Hj}]} else {lappend retlist 0.}
       
        # Output the types of correlated nucleotides
        #vputs [format "%[expr {$lenListNt*7-7}]s Pairs =" " "] 7 1
        for {set i 0} {$i<=$numSignificantPair} {incr i} {
            #vputs " [lindex $significantPair($i) 1]:[lindex $significantPair($i) 2]" 7 1
            lappend pairs [lindex $significantPair($i) 1]:[lindex $significantPair($i) 2]
        }
        #vputs " " 7
        lappend retlist $pairs
    
        return $retlist
    }
    # upe


      
    ###   SetupChi2Distr   
    #
    # Table of chi-square probabilities
    #
    #	For chi^2 analysis see for example
    #		Duerr, W. & Mayer, H. (1992).
    #		Wahrscheinlichkeitsrechnung und schließende Statistik.
    #		Carl Hanser Verlag, M"unchen.
    #
    #		Anderson, T.W. & Sclove, S.L. (1978).
    #		An Introduction to the Statistical Analysis of Data.
    #		Houghton Mifflin Company, Boston.
    #
    proc SetupChi2Distr {} {
    ###################
   
       # Approximate quantiles of the chi-square distribution
       # 
       #      |                                  p
       #   df |   .005    .01     .025   .05      .10    .90      .95     .975    .99     .995
       #  ----+-------------------------------------------------------------------------------
       #    1 |   .00004  .00016  .00098  .0039   .0158  2.71    3.84    5.02    6.63    7.88
       #    2 |   .0100   .0201   .0506   .1026   .2107  4.61    5.99    7.38    9.21   10.60
       #    3 |   .0717   .115	   .216	   .352    .584   6.25    7.81    9.35   11.34   12.84
       #    4 |   .207	   .297	   .484	   .711   1.064   7.78    9.49   11.14   13.28   14.86
       #    5 |   .412	   .554	   .831   1.15    1.61    9.24   11.07   12.83   15.09   16.75
       #    6 |   .676    .872   1.24    1.64    2.20   10.64   12.59   14.45   16.81   18.55
       #    7 |   .989   1.24    1.69    2.17    2.83   12.02   14.07   16.01   18.48   20.28
       #    8 |  1.34    1.65    2.18    2.73    3.49   13.36   15.51   17.53   20.09   21.96
       #    9 |  1.73    2.09    2.70    3.33    4.17   14.68   16.92   19.02   21.67   23.59
       #   10 |  2.16    2.56    3.25    3.94    4.87   15.99   18.31   20.48   23.21   25.19
       #   11 |  2.60    3.05    3.82    4.57    5.58   17.28   19.68   21.92   24.73   26.76
       #   12 |  3.07    3.57    4.40    5.23    6.30   18.55   21.03   23.34   26.22   28.30
       #   13 |  3.57    4.11    5.01    5.89    7.04   19.81   22.36   24.74   27.69   29.82
       #   14 |  4.07    4.66    5.63    6.57    7.79   21.06   23.68   26.12   29.14   31.32
       #   15 |  4.6     5.23    6.26    7.26    8.55   22.31   25	 27.49   30.58   32.80
       #   16 |  5.14    5.81    6.91    7.96    9.31   23.54   26.30   28.85   32.00   34.27
       #   18 |  6.26	  7.01	  8.23    9.39   10.86   25.99   28.87   31.53   34.81   37.16
       #   20 |  7.43	  8.26	  9.59   10.85   12.44   28.41   31.41   34.17   37.57   40.00
       #   24 |  9.89	 10.86   12.40   13.85   15.66   33.20   36.42   39.36   42.98   45.56
       #   30 | 13.79   14.95   16.79   18.49   20.60   40.26   43.77   46.98   50.89   53.67
       #   40 | 20.71   22.16   24.43   26.51   29.05   51.81   55.76   59.34   63.69   66.77
       #   60 | 35.53   37.48   40.48   43.19   46.46   74.40   79.08   83.30   88.38   91.95
       #  120 | 83.85   86.92   91.58   95.70  100.62  140.23  146.57  152.21  158.95	163.64
       
       #    set chi_values [list \
       #   .00004  .00016  .00098  .0039   .0158  2.71    3.84    5.02    6.63    7.88 \
       #   .0100   .0201   .0506   .1026   .2107  4.61    5.99    7.38    9.21   10.60 \
       #   .0717   .115    .216    .352    .584   6.25    7.81    9.35   11.34   12.84 \
       #   .207    .297    .484    .711   1.064   7.78    9.49   11.14   13.28   14.86 \
       #   .412    .554    .831   1.15    1.61    9.24   11.07   12.83   15.09   16.75 \
       #   .676    .872   1.24    1.64    2.20   10.64   12.59   14.45   16.81   18.55 \
       #   .989   1.24    1.69    2.17    2.83   12.02   14.07   16.01   18.48   20.28 \
       #  1.34    1.65    2.18    2.73    3.49   13.36   15.51   17.53   20.09   21.96 \
       #  1.73    2.09    2.70    3.33    4.17   14.68   16.92   19.02   21.67   23.59 \
       #  2.16    2.56    3.25    3.94    4.87   15.99   18.31   20.48   23.21   25.19 \
       #  2.60    3.05    3.82    4.57    5.58   17.28   19.68   21.92   24.73   26.76 \
       #  3.07    3.57    4.40    5.23    6.30   18.55   21.03   23.34   26.22   28.30 \
       #  3.57    4.11    5.01    5.89    7.04   19.81   22.36   24.74   27.69   29.82 \
       #  4.07    4.66    5.63    6.57    7.79   21.06   23.68   26.12   29.14   31.32 \
       #  4.6     5.23    6.26    7.26    8.55   22.31   25      27.49   30.58   32.80 \
       #  5.14    5.81    6.91    7.96    9.31   23.54   26.30   28.85   32.00   34.27 \
       #  6.26    7.01    8.23    9.39   10.86   25.99   28.87   31.53   34.81   37.16 \
       #  7.43    8.26    9.59   10.85   12.44   28.41   31.41   34.17   37.57   40.00 \
       #  9.89   10.86   12.40   13.85   15.66   33.20   36.42   39.36   42.98   45.56 \
       # 13.79   14.95   16.79   18.49   20.60   40.26   43.77   46.98   50.89   53.67 \
       # 20.71   22.16   24.43   26.51   29.05   51.81   55.76   59.34   63.69   66.77 \
       # 35.53   37.48   40.48   43.19   46.46   74.40   79.08   83.30   88.38   91.95 \
       # 83.85   86.92   91.58   95.70  100.62  140.23  146.57  152.21  158.95  163.64]
       # 
       #    set chi2_p [list 0.005 0.010 0.025 0.050 0.100 0.900 0.950 0.975 0.990 0.995]
       #    set index 0
       #    for {set i 1} {$i<24} {incr i} {;	# degrees of freedom
       #       if       {$i==17} {
       #          set df 18
       #       } elseif {$i==18} {
       #          set df 20
       #       } elseif {$i==19} {
       #          set df 24
       #       } elseif {$i==20} {
       #          set df 30
       #       } elseif {$i==21} {
       #          set df 40
       #       } elseif {$i==22} {
       #          set df 60
       #       } elseif {$i==23} {
       #          set df 120
       #       } else {
       #          set df $i
       #       }
       #       foreach p $chi2_p {;	# probability
       #          set chi2($df,$p) [lindex $chi_values $index]
       # 	 incr index
       #       }
       #    }
       #    foreach p $chi2_p {
       #       set chi2(17,$p) $chi2(18,$p)
       #       set chi2(19,$p) $chi2(20,$p)
       #       for {set i 21} {$i<24} {incr i} {
       #          set chi2($i,$p) $chi2(24,$p)
       #       }
       #       for {set i 25} {$i<30} {incr i} {
       #          set chi2($i,$p) $chi2(30,$p)
       #       }
       #       for {set i 31} {$i<40} {incr i} {
       #          set chi2($i,$p) $chi2(40,$p)
       #       }
       #       for {set i 41} {$i<60} {incr i} {
       #          set chi2($i,$p) $chi2(60,$p)
       #       }
       #       for {set i 61} {$i<120} {incr i} {
       #          set chi2($i,$p) $chi2(120,$p)
       #       }
       #    }
       #    set chi2_p [list .995 .990 .975 .950 .900 .100 .050 .025 .010 .005]
       
       # Approximate quantiles of the chi-square distribution (Bronstein)
       # 
       #      |                          p
       #   df |   .01   .02  .05    .10   .20   .30   .50   .70   .80   .90   .95   .98   .99   .995   .998   .999
       #  ----+---------------------------------------------------------------------------------------------------
       #    9 |  2.09  2.53  3.32  4.17  5.38  6.39  8.34 10.7  12.2  14.7  16.9  19.7  21.7  23.6   26.1   27.9
       #   16 |  5.8   6.6   8.0   9.3  11.2  12.6  15.3  18.4  20.5  23.5  26.3  29.6  32.0  34.3   37.1   39.3
       #   25 | 11.5  12.7  14.6  16.5  18.9  20.9  24.3  28.2  30.7  34.4  37.7  41.6  44.3  46.9   50.1   52.6
       
       	variable chi2
       	variable chi2_p
    		
       	variable debug_level
    		variable verbose_level
    		
    		
       	if {$debug_level>1} {
            puts "[namespace current]::SetupChi2Distr executing"
        } elseif {$verbose_level>1} {
            puts "Setting up Chi^2 Distribution"
        }
       	
       	set chi2_p [list .01 .02 .05 .10 .20 .30 .50 .70 .80 .90 .95 .98 .99 .995 .998 .999]
       	set chi_values [list \
       	       2.09  2.53  3.32  4.17  5.38  6.39  8.34 10.7  12.2  14.7  16.9  19.7  21.7  23.6   26.1   27.9\
       	       5.8   6.6   8.0   9.3  11.2  12.6  15.3  18.4  20.5  23.5  26.3  29.6  32.0  34.3   37.1   39.3\
       	      11.5  12.7  14.6  16.5  18.9  20.9  24.3  28.2  30.7  34.4  37.7  41.6  44.3  46.9   50.1   52.6\
       	    ]
       	set index 0
       	foreach df [list 9 16 25] {     ; # degrees of freedom
       		foreach p $chi2_p {         ; # probability
       			set chi2($df,$p) [lindex $chi_values $index]
       	 		incr index
       	   }
       	}
       	set chi2_p [list .999 .998 .995 .99 .98 .95 .90 .80 .70 .50 .30 .20 .10 .05 .02 .01]
   }
   # SetupChi2Distr
   
   
   
   
    ###   UnsetChi2Distr   
    #
    #     Unsets the Chi2 Distribution to save some memory
    #     Normally not needed
    #
    # Arguments:
    #     None
    #
    # Results:
    #     Namespace wide chi2 distribution is unset
    #
    proc UnsetChi2Distr {} {
    ###################
   	  variable chi2
   	  variable chi2_p
   	
        catch {unset chi2_p}
        catch {unset chi2}
    }
    # UnsetChi2Distr




    
    ###   CheckPlausibility   
    #
    #     Checks plausibility of values
    #
    #
    # Arguments:
    #     !!! TODO !!!
    #
    # Results:
    #     Returns ERROR if inconsistenty is detected and prints out some info
    #     OK otherwise
    #
    #
    proc CheckPlausibility {listNt i j n ff} {
    #######################
    	  variable debug_level
    	  variable verbose_level
    	
    	upvar $i fi
    	upvar $j fj
    	upvar $ff f
    	
    	set sumi 0.
        set sumj 0.
    	set total 0.
    	
    	
   	
        foreach a $listNt {
            if {$fi($a)<0. || $fi($a)>$n} {
                puts "[namespace current]::upe: Illegal size of a frequency(i)!"
                return ERROR
            }
            if {$fj($a)<0. || $fj($a)>$n} {
                puts "[namespace current]::upe: Illegal size of a frequency(j)!"
                return ERROR
   		    }
   		    set sumi [expr {$sumi + $fi($a)}]
   		    set sumj [expr {$sumj + $fj($a)}]
   	    }
    	
   	    if {$sumi!=$n || $sumj!=$n} {
   		    puts "[namespace current]::upe: The sum of single frequencies is inconsistent!"
   		    puts "[namespace current]::upe:    sum(i) = $sumi; sum(j) = $sumj; n = $n"
   		    return ERROR
   	    }
   	
   	
   	    foreach a $listNt {
   		    set sum 0.
   		    foreach b $listNt {
   			    if {$f($a,$b)<0. || $f($a,$b)>$fi($a)} {
   				    puts "[namespace current]::upe: Illegal size of a frequency(I,j)!"
   				    return ERROR
   			    }
   			    set sum   [expr {$sum   + $f($a,$b)}]
   			    set total [expr {$total + $f($a,$b)}]
   		    }
  			if {$sum!=$fi($a)} {
  				puts "[namespace current]::upe: The sum of frequencies(I,j) is inconsistent!"
  				return ERROR
  			}
  		}
    	
  		if {$total!=$n} {
  			puts "[namespace current]::upe: The total sum of frequencies(i,j) is inconsistent!"
  			return ERROR
  		}
   
    
  		foreach b $listNt {
  			set sum 0.
  			foreach a $listNt {
  				if {$f($a,$b)<0. || $f($a,$b)>$fj($b)} {
  					puts "[namespace current]::upe: Illegal size of a frequency(i,J)!" 7
  					return ERROR
  				}
  				set sum [expr {$sum + $f($a,$b)}]
  			}
  			if {$sum!=$fj($b)} {
  				puts "[namespace current]::upe: The sum of frequencies(i,J) is inconsistent!" 7
  				return ERROR
  			}
  		}
    	
    	return OK
    }
    # CheckPlausibility
    
    
    
    
    ###   SetVar   
    #
    #     Use to set namespace wide variables
    #     Alternativ to "namespace eval set variable value"
    #
    # Arguments:
    #     varname:: name of namespace wide variable
    #     val     : value of the variable
    #     Syntax is basically the same as "set"
    #
    # Results
    #     Namespace variable value is changed    
    #
    proc SetVar {varname val} {
    	variable $varname
    	set $varname $val	
    }
    # SetVar
    
};# namespace eval 
