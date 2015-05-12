namespace eval logo {

	variable logoInput
	variable aInput
	variable cInput
	variable gInput
	variable uInput
	variable auInput
	variable cgInput
	variable guInput
	variable startPos
	variable struct
	variable field
	variable type
	variable useZero
	variable useT
	variable struct_alnnum
		variable workdir
	global seq

	# ---   updateLogo
	#
	#
	#
	#Update logoInput
	proc updateLogo {} {
	
		variable logoInput
		variable aInput
		variable cInput
		variable gInput
		variable uInput
		variable auInput
		variable cgInput
		variable guInput
		variable startPos
		variable struct
		variable field
		variable type
		variable useZero
		variable useT
		variable struct_alnnum
	 	variable title
		
		global seq
	
#        puts "\nntprobs:       $aInput++$cInput++$gInput++$uInput"
#        puts "probAU:        $auInput"
#        puts "probCG:        $cgInput"
#        puts "probGU:        $guInput"
#        puts "strucseq:      $struct"
#        puts "assignment:    $field"
#        puts "logotype:      $type"
#        puts "startpos:      $startPos"
#        puts "poslabel:      $useZero"
#        puts "TorUs:         $useT"
#        puts "struct_alnnum:\n$struct_alnnum"
#	    for {set s 1} {$s<=$seq(n_seq)} {incr s} {
#	        puts [string toupper $seq(nt,$s)]
#	    }

		set logoInput "ntprobs=$aInput++$cInput++$gInput++$uInput&probAU=$auInput&probCG=$cgInput&probGU=$guInput&strucseq=$struct&assignment=$field&logotype=$type&startpos=$startPos&poslabel=$useZero&TorUs=$useT&comments="
		#puts $logoInput
		append logoInput "%3E+";    # "> "
	   #append logoInput [StructTransform $bplist $seq(aln_len)]
		append logoInput $struct_alnnum
	   append logoInput "%0D%0A";  # "<cr>"
	   for {set s 1} {$s<=$seq(n_seq)} {incr s} {
	        append logoInput "%3E+"
	        append logoInput [string toupper $seq(nt,$s)]
	        append logoInput "%0D%0A"
	    }
	
		logo::upload 
	
	}



	# --- frontEnd
	#
	#
	#
	proc frontEnd {struct_arg name dir} {
		#
		variable logoInput
		variable aInput
		variable cInput
		variable gInput
		variable uInput
		variable auInput
		variable cgInput
		variable guInput
		variable startPos
		variable struct
		variable field
		variable type
		variable useZero
		variable useT
		variable struct_alnnum
		variable title
		variable workdir

        
		set title $name		
		set struct_alnnum $struct_arg
		set workdir $dir

        set win .structlogowin
        catch {destroy $win}
        toplevel    $win
        wm title    $win "Structure Logo"
        wm iconname $win StructureLogo
        
        
        
		#Buttons zum Beenden uns Abschicken
		set buttonFrame [frame $win.buttons] 
		set gobutton [button $win.buttons.gobutton -text "Send" -command "logo::updateLogo;destroy $win"]
		set quitbutton [button $win.buttons.quitbutton -text "Quit" -command "destroy $win"]
        set infotxt "(By clicking \"Send\" you request a Logo from http://www.cbs.dtu.dk/~gorodkin/appl/slogo.html)"
        set infolbl [label $win.buttons.message -justify center -text "$infotxt"]
        #
        grid $win.buttons.gobutton -row 0 -column 0
        grid $win.buttons.quitbutton -row 0 -column 1
        grid $win.buttons.message -row 1 -column 0 -columnspan 2

		

		#Eingabe der NT-Wahrscheinlichkeit; Summe muss immer 1.00 sein
		set ntInputFrame [frame $win.ntInput] 
		frame $win.ntInput.aFrame 
		frame $win.ntInput.cFrame 
		frame $win.ntInput.gFrame 
		frame $win.ntInput.uFrame 
		
		label $win.ntInput.title -text "NT Prob:"
		set aNt [entry $win.ntInput.aFrame.aNt -textvariable logo::aInput -width 6]
		set aPrompt [label $win.ntInput.aFrame.aPrompt -text "A:"]
		set cNt [entry $win.ntInput.cFrame.cNt -textvariable logo::cInput -width 6]
		set cPrompt [label $win.ntInput.cFrame.cPrompt -text "C:"]
		set gNt [entry $win.ntInput.gFrame.gNt -textvariable logo::gInput -width 6]
		set gPrompt [label $win.ntInput.gFrame.gPrompt -text "G:"]
		set uNt [entry $win.ntInput.uFrame.uNt -textvariable logo::uInput -width 6]
		set uPrompt [label $win.ntInput.uFrame.uPrompt -text "U:"]
		
		set aInput "0.25"
		set cInput "0.25"
		set gInput "0.25"
		set uInput "0.25"
		
		pack $win.ntInput.title -side top 
		pack $win.ntInput.aFrame -side top
		pack $win.ntInput.cFrame -side top
		pack $win.ntInput.gFrame -side top
		pack $win.ntInput.uFrame -side top
		
		pack $win.ntInput.aFrame.aNt -side right
		pack $win.ntInput.aFrame.aPrompt -side left
		pack $win.ntInput.cFrame.cNt -side right
		pack $win.ntInput.cFrame.cPrompt -side left
		pack $win.ntInput.gFrame.gNt -side right
		pack $win.ntInput.gFrame.gPrompt -side left
		pack $win.ntInput.uFrame.uNt -side right
		pack $win.ntInput.uFrame.uPrompt -side left
		
		#
		#
		#Basenpaargewichtung
		set bpInputFrame [frame $win.bpInput]
		set auFrame [frame $win.bpInput.auFrame]
		set cgFrame [frame $win.bpInput.cgFrame]
		set guFrame [frame $win.bpInput.guFrame]
		
		label $win.bpInput.title -text "BP weight:"
		set auW [entry $win.bpInput.auFrame.auW -textvariable logo::auInput -width 6]
		set auPrompt [label $win.bpInput.auFrame.auPrompt -text "AU:"]
		set cgW [entry $win.bpInput.cgFrame.cgW -textvariable logo::cgInput -width 6]
		set cgPrompt [label $win.bpInput.cgFrame.cgPrompt -text "CG:"]
		set guW [entry $win.bpInput.guFrame.guW -textvariable logo::guInput -width 6]
		set guPrompt [label $win.bpInput.guFrame.guPrompt -text "GU:"]
		
		set auInput "1.0"
		set cgInput "1.0"
		set guInput "1.0"
		
		pack $win.bpInput.title -side top
		pack $win.bpInput.auFrame -side top
		pack $win.bpInput.cgFrame -side top
		pack $win.bpInput.guFrame -side top
		
		pack $win.bpInput.auFrame.auW -side right
		pack $win.bpInput.auFrame.auPrompt -side left
		pack $win.bpInput.cgFrame.cgW -side right
		pack $win.bpInput.cgFrame.cgPrompt -side left
		pack $win.bpInput.guFrame.guW -side right
		pack $win.bpInput.guFrame.guPrompt -side left
		
		#
		#
		#Start Position
		set  ntStartPos [frame $win.ntStart]
		
		entry $win.ntStart.ntPos -textvariable logo::startPos -width 6
		label $win.ntStart.ntPosLbl -text "Start Position:"
		set startPos "1"
		grid $win.ntStart.ntPosLbl -row 0
		grid $win.ntStart.ntPos -row 1
		
		#
		#
		#StructureLogo
		set position 0
		set row 1
		set logoList [list "Structure" Y "Sequence" N]
		set structFrame [frame $win.struct]
		set struct Y
        
        label $win.struct.lbl -text "Logo Type"
        grid $win.struct.lbl -row 0 -column 0
        incr row
		foreach {item logo} $logoList {
			radiobutton $win.struct.$position  \
				-text "$item"  -variable logo::struct -value $logo   
            grid $win.struct.$position -column 0 -row $row
            incr position
            incr row
        }
		#
		#
		#Logo Type
		set position 0
		set row 1
		set logoType [list "Logo Type 1" 1 "Logo Type 2" 2]
		set logoFrame [frame $win.logo]
		set type 2
		
		foreach {item lType} $logoType {
		
			radiobutton $win.logo.$position \
				-text "$item" -variable logo::type -value $lType
				grid $win.logo.$position -column 1 -row $row
				incr position
				incr row
			}

        
		#
		#
		#Field Assignment
		set position 0
		set row 1
		set assignmentList [list "Assignment Field" Y "No Assignment Field" N]
		set fieldFrame [frame $win.field]
		set field Y
		
		foreach {item assignment} $assignmentList {
		
			radiobutton $win.field.$position \
				-text "$item" -variable logo::field -value $assignment
				grid $win.field.$position -column 1 -row $row -sticky w
				incr position
				incr row
		
			}
		
		
		#
		#
		#Use Zero...
		set position 0
		set row 1
		set zeroList [list "Use 0 as start" Y "Use 1 as start" N]
		set zeroFrame [frame $win.zero]
		set useZero Y
		
		foreach {item useZ} $zeroList {
		
			radiobutton $win.zero.$position \
				-text "$item" -variable logo::useZero -value $useZ
				grid $win.zero.$position -column 1 -row $row -sticky w
				incr position
				incr row
		
			}
		
		#
		#
		#Display "T"
		set position 0
		set row 1
		set tList [list "Display U" N "Display T" Y]
		set tFrame [frame $win.t]
		set useT N
		
		foreach {item t} $tList {
		
			radiobutton $win.t.$position \
				-text "$item" -variable logo::useT -value $t
				grid $win.t.$position -column 1 -row $row
				incr position
				incr row
		
			}
		
		#
		#
		#Erstellen des Widgets
# 		pack $ntInputFrame -side top -side left
# 		pack $bpInputFrame -side top -side right
# 		pack $ntStartPos -side top  
# 		pack $win.buttons -side bottom
# 		pack $structFrame  -side left
# 		pack $zeroFrame -side left 
# 		pack $tFrame -side left 
# 		pack $fieldFrame  -side left
# 		pack $logoFrame  -side left
 		grid $ntInputFrame -row 0 -column 0 -sticky n
 		grid $bpInputFrame -row 0 -column 1 -sticky n
 		grid $ntStartPos -row 0 -column 2 -sticky n
 		grid $structFrame -row 0 -column 3 -sticky n
 		grid $zeroFrame -row 0 -column 4 -sticky n
 		grid $tFrame -row 0 -column 5 -sticky n
 		grid $fieldFrame -row 0 -column 6 -sticky n
 		grid $logoFrame -row 0 -column 7 -sticky n
 		grid $win.buttons -row 1 -column 0 -columnspan 8
        

        raise $win
        focus $win

	}


	# --- upload
	#
	#
	#
	proc upload {} {
	
		variable logoInput
		variable aInput
		variable cInput
		variable gInput
		variable uInput
		variable auInput
		variable cgInput
		variable guInput
		variable startPos
		variable struct
		variable field
		variable type
		variable useZero
		variable useT
		variable struct_alnnum
		variable title
		variable workdir
	   global seq

		#set struct_alnnum $struct_arg
		#set title $name

		# Alles, was jetzt hier folgt, sollte in eine extra Prozedur 
		#   gepackt werden, da das nichts mehr mit ConStruct zu
		#   tun hat
	    set unique [clock clicks]
	    set outFN      [file join / tmp "${title}_${unique}.out"];       # Hierin sollte die Antwortseite stehen
	    set stderrbuff [file join / tmp "${title}_${unique}.stderr"];    # Instead of writing to stderr write to file
	    if {[file exists $outFN]}      {file delete $outFN}
	    if {[file exists $stderrbuff]} {file delete $stderrbuff}
	    set return [ \
	        curl::transfer -url http://www.cbs.dtu.dk/cgi-bin/rnalogoform   \
	            -verbose        1           \
	            -post           1           \
	            -postfields     $logoInput  \
	            -file           $outFN      \
	            -errorbuffer    errorbuff   \
	            -stderr         $stderrbuff \
	    ]
	    if {$return!=0} {
	        puts "ERROR: Can't produce logo; error code $return"
	        puts "$errorbuff"
	        puts "       More information might be in $stderrbuff"
	        return
	    } else {
	        if {[file exists $stderrbuff]} {file delete $stderrbuff}
	    }
	
	    # put together the URL for PostScript file
	    set ps_url "dummy"	
	    set outFH [open $outFN r]
	        while {[gets $outFH line]>=0} {
	            if {[regexp "<A HREF=\"(.*)html\"" $line schrott ps_url]} {
	    	        regsub "Sub" $ps_url "" ps_url
	    	        append ps_url "ps"
	    	        break
	            }
	        }
	    close $outFH
	    file delete $outFN
	    if {$ps_url=="dummy"} {
	        puts "ERROR: Can't get URL of postscript file"
	        return
	    }

        set outFN [tk_getSaveFile -title "Select structure file" \
		  		-initialdir $workdir -initialfile "${title}_logo.ps"]
        if {$outFN==""} {
				puts "Save dialog aborted..."
				return
        }

	
	    # get PostScript file
	    #set outFN [file join $workdir "${title}_logo.ps"];  # PostScript-Datei
	    if {[file exists $outFN]}       {file delete $outFN}
	    if {[file exists $stderrbuff]}  {file delete $stderrbuff}
	    set return [ \
	        curl::transfer -url $ps_url     \
	            -file           $outFN      \
	            -errorbuffer    errorbuff   \
	            -stderr         $stderrbuff \
	    ]
	    if {$return!=0} {
	        puts "ERROR: Can't get PostScript file; error code $return"
	        puts "$errorbuff"
	        puts "       More information might be in $stderrbuff"
	        return
	    } else {
	        if {[file exists $stderrbuff]} {file delete $stderrbuff}
	    }
	
	    # Test for a PostScript file
	    #   otherwise there was an error
	    set doit "eval exec file -b $outFN"
	    set ps "PostScript"
	    if {[catch $doit fehler]==0 && \
	        [string compare -length [string length $ps] $fehler $ps]!=0} {
	        puts "ERROR: Got wrong file type ($fehler)"
	        return
	    }
	
	    # view PostScript file
	    #   assume that the command ends with '-'
	    global opts
	    set dummy "${opts(print_cmd,screen)}-orientation=landscape"
	    catch {eval exec $dummy $outFN &} result
	}
	
	
}; # namespace eval logo

