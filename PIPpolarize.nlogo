extensions [ vid ]

globals [
  setup-success?
  run-index

  ; Internally set Constants ;
  colorname-outpatches
  RGB-pip1
  RGB-pip2
  RGB-kinase
  RGB-pptase

  ; Calculated Constants ;
  dist_enz ; distance travelled by enzyme per timestep
  alpha-pip ; alpha = beta (from FTCS calculation)
  alpha-enz ; alpha = beta (from FTCS calculation)
  patchLength
  const-k_Poff
  const-p_Poff

  Pstay

  ; Patch lists ;
  inpatches
  outpatches
  Lpatches
  Spatches

  ; Trackers - ever changing ;
  avg_x
  time
  next_tlapse_time

  file-prefix
]

breed [ kinases kinase ]
breed [ pptases pptase ]
breed [ clocks clock ]

patches-own [
  real_neighbors
  n_neighbors
  x_patch
  dx_patch
  updated_xpatch ; this is a must-have

  ; Deterministic mode variables
  patch_k_density
  patch_p_density
  updated_kpatch
  updated_ppatch

  ; Stochastic mode variables
  k_Pon
  p_Pon
]

to setup
  clear-all
  display
  set setup-success? true ; It will change to false if something goes wrong in the resetting process.

  RESET-TICKS

  plot_dxdt_vs_x

  ; Fix constants (internally defined) ;
  set colorname-outpatches brown - 1
  set RGB-pip1 [0 100 255] ;  set RGB-pip1 extract-rgb blue
  set RGB-pip2 [255 200 0]
  set RGB-kinase [255 150 0]
  set RGB-pptase [0 150 255]

  ; Init patchsets
  set inpatches no-patches
  set outpatches no-patches

  ; Setup time
  set time 0

  ; Setup space
  resize-world 0 (nGrid - 1) 0 (nGrid - 1)
  if smaller-part-test? and Calculation-Type = "deterministic" and geometry-setup = "None" [resize-world 0 (nGrid - 1) 0 (nGrid - 1) / 10]
  set-patch-size world_pixel_length / nGrid
  let wrap? true
  ifelse wrap? [ __change-topology true  true  ]
                [ __change-topology false false ]
  setup_world_from_input_file
  set patchLength worldLength / nGrid
  set-neighbors_and_outColors
  initialize_patches

  ; L patch and S patch
  if plot-xL-xS? [
    set Lpatches no-patches
    set Spatches no-patches
    ask inpatches [
      if pxcor >= 24 / 50 * nGrid and pxcor < 26 / 50 * nGrid and pycor >= 24 / 50 * nGrid and pycor < 26 / 50 * nGrid [ set Lpatches (patch-set Lpatches self)]
      if pxcor >= 45 / 50 * nGrid and pxcor < 47 / 50 * nGrid and pycor >= 44 / 50 * nGrid and pycor < 46 / 50 * nGrid [ set Spatches (patch-set Spatches self)]
    ]
    if Spatches = no-patches or Lpatches = no-patches [ user-message "Are you using 50-snail6.png as geometry? If not, try setting the \"plot-xL-xS? switch\" to \"off\"." ]
  ]

  ; Set-up clock
  if timestamp-on-image?
  [  create-clocks 1 [ set label (precision time 1) setxy (min-pxcor + nGrid / 5) (max-pycor - nGrid / 20) set size 0 ]  ]

  ; Setup (unchanging) constants
  set alpha-pip D_pip * timestep / (worldLength / nGrid) ^ 2
  if Calculation-Type = "deterministic" [set alpha-enz D_enz * timestep / (worldLength / nGrid) ^ 2]
  set Pstay 1 - (4 * D_enz * timestep) / (worldLength / nGrid) ^ 2
  set const-k_Poff k_koff * timestep   if const-k_Poff > 1 [user-message "k_off_prob_sto > 1"]
  set const-p_Poff p_koff * timestep   if const-p_Poff > 1 [user-message "p_off_prob_sto > 1"]

  ; Check for stability condition for the numerical treatment of diffusion
  if timestep > (worldLength / nGrid) ^ 2 / (4 * D_pip) [ user-message (word "FTCS dispersion - Stable condition not met. Decrease timestep below " ((worldLength / nGrid) ^ 2 / (4 * D_pip)) ) ]
  if Calculation-Type = "deterministic" [if timestep > (worldLength / nGrid) ^ 2 / (4 * D_enz) [ user-message (word "Stable condition not met. Decrease timestep below " ((worldLength / nGrid) ^ 2 / (4 * D_enz)) ) ]]

  ; Init save files
  set file-prefix retrieve-simul_info_string
  ifelse save_timelapse_img? or record_vid? or save_all_plots? or save-xL-xS?  [  set setup-success? initialize-saving   ]
  [ set save-dir-name "N/A" ]
  display
end


to plot_dxdt_vs_x
  let Nxpoints 100
  let xpoints n-values (Nxpoints + 1) [ x -> x / 100 ]
  set-current-plot "dxdt vs x"

  carefully [
    if Enzyme-Pair-Type = "memK-solP" [
      foreach xpoints [
        x -> let dxdt 0
        set dxdt k_mkon / k_koff * x * k_mkcat * (1 - x) / (k_mKm + 1 - x) - solP_mkcat * x / (p_mKm + x)
        plotxy x dxdt
      ]
      set-current-plot-pen "dxdt0"
      foreach xpoints [x -> plotxy x 0]
    ]
    if Enzyme-Pair-Type = "memK-memP" [
      foreach xpoints [
        x -> let dxdt 0
        set dxdt k_mkon / k_koff * x * k_mkcat * (1 - x) / (k_mKm + 1 - x) - p_mkon / p_koff * (1 - x) * memP_mkcat * x / (p_mKm + x)
        plotxy x dxdt
      ]
      set-current-plot-pen "dxdt0"
      foreach xpoints [x -> plotxy x 0]
    ]
    set-current-plot-pen "x0.5"
    plotxy 0.5 plot-y-max
    plotxy 0.5 plot-y-min
  ] [
    print (word "dx/dt vs x Plotting failed (check kinetic parameters) Code:" error-message)
  ]
end


to initialize_patches
  ; Read input k-, p- density and x.
  let init-k_patch-density item 0 read-from-string KIN-PPT-X
  let init-p_patch-density item 1 read-from-string KIN-PPT-X
  if Enzyme-Pair-Type = "memK-solP" and init-p_patch-density != 0
  [ user-message "**Warning: \"Setting up a solution-phosphatase setup\" and still adding binding phosphatase on the membrane as initial condition?"]
  let init-x_patch item 2 read-from-string KIN-PPT-X

  ; Clear out enzymes on the membrane.
  ask kinases [die]
  ask pptases [die]
  ask inpatches [
    set patch_k_density 0
    set patch_p_density 0
  ]

  ; Placing kinases and phosphatases
  if Calculation-Type = "deterministic" [
    ask inpatches [
      set patch_k_density init-k_patch-density
      set patch_p_density init-p_patch-density
    ]
  ]

  if Calculation-Type = "stochastic" [
    let n_init-kinases init-k_patch-density * count inpatches * patchLength ^ 2
    repeat n_init-kinases[
      ask one-of inpatches [
        sprout-kinases 1
        [
          set shape "circle"
          set size enz_size
          set color RGB-kinase
    ]  ]  ]
    let n_init-pptases init-p_patch-density * count inpatches * patchLength ^ 2
    repeat n_init-pptases[
      ask one-of inpatches [
        sprout-pptases 1
        [
          set shape "circle"
          set size enz_size
          set color RGB-kinase
    ]  ]  ]
  ]

  ; Setting up x (equilvalent to setting up PIP1/PIP2 ratio)
  ask inpatches [   set x_patch init-x_patch   ]
  if Calculation-Type = "deterministic" [ fluctuate ] ; If kinetics is deterministic, initial fluctuation is "required" to see any kind of instability to develop

  ; Visualize x as patch color
  ask inpatches [ represent-x-as-patch-color ]

  set avg_x mean [x_patch] of inpatches

  ; Show-or-hide enzymes
  if Calculation-Type = "stochastic" [
    if show_enz? [ ask kinases [show-turtle] ask pptases [show-turtle] ]
    if not show_enz? [ ask kinases [hide-turtle] ask pptases [hide-turtle] ]
  ]
end

to-report initialize-saving
  let success? true ; This will change to false if something goes wrong in the process.
  if (save_timelapse_img? or record_vid? or save_all_plots? or save-xL-xS?) and (save-dir-name = "N/A" or save-dir-name = "") [
    print "wft"
    carefully [      set save-dir-name user-directory    ]
    [ ;user-message "Results saving directory not selected."
      set success? false ; This will fall-into the "Saving to the save-dir-name failed." case in the down below.
      set save-dir-name "N/A"
    ]
  ]
  if record_vid? [ carefully [vid:start-recorder]
    [ user-message "Video init failed. Select a directory to dump the video (if possible)."
      set-current-directory user-directory
      carefully [vid:save-recording (word "dump_mov.mp4")]
      [ user-message "If there was any content, the video is dumped as dump_mov.mp4."]
      set success? false
    ]
  ]
  carefully [ export-interface (word save-dir-name "iface-t0 " file-prefix ".png") ] ; If this directory does not exist, this will spit out an error message.
    [ user-message "Saving to the save-dir-name failed. Make sure you've put in a valid directory"
      set success? false
      set save-dir-name "N/A"
  ]
  report success?
end


to go
  if setup-success? = false
  [ user-message "Reset status unsuccessful."    stop  ]
  if save-dir-name != "N/A" and not file-exists? (word save-dir-name "iface-t0 " file-prefix ".png") ; Trying to find the screenshot of the initial interface taken at "setup".
  [ user-message "save-dir-name (probably) changed since \"setup\"" set save-dir-name "N/A" set setup-success? false stop ] ; If it's not found, tell the user that probably, you changed the save-dir after you pressed "setup"


  if "no" = user-one-of (word "Saving result file(s) to " save-dir-name " - Okay?") ["yes" "no"]
    [stop]



  ; Check stop-conditions and apply necessary ending steps
  if time > endtime
  [
    ; Save
    if save_timelapse_img? or record_vid? or save_all_plots? or save-xL-xS?
    [
      set-current-directory save-dir-name
      export-interface (word "iface-End " file-prefix ".png")
    ]

    set run-index  run-index + 1

    ; single-run ends and moves on to the next run
    if run-index < N-runs [
      ask kinases [ die ]
      ask pptases [ die ]
      initialize_patches
      set time 0
      set next_tlapse_time 0
    ]

    ; All-runs end and Export "xL-xS"
    if run-index >= N-runs [
      if save-xL-xS? [export-plot "xL-xS" (word "xL-xS of " file-prefix  " " N-runs "-runs.csv")]
      if record_vid? [  vid:save-recording (word file-prefix "_mov.mp4") ]
      if save_all_plots? [export-all-plots (word file-prefix " - allplots.csv")]
      set run-index  run-index - 1 ; for visual purpose
      stop
    ]
    ; Clear non-accumulative plots
    set-current-plot "Max Min dx_patch" clear-plot
    set-current-plot "Max Pon-patch" clear-plot
;    set-current-plot "Max - Min dx_patch" clear-plot
  ]

  ; Visual/graph updates
  if show_enz? [ ask kinases [show-turtle] ask pptases [show-turtle] ]
  if not show_enz? [ ask kinases [hide-turtle] ask pptases [hide-turtle] ]
  if timestamp-on-image? [ ask clocks [ set label (precision time 1) ] ]
  ask inpatches [ represent-x-as-patch-color ]
  set avg_x mean [x_patch] of inpatches
  if plot-xL-xS? [
    set-current-plot "xL-xS"
    set-current-plot-pen (word (run-index))
    plotxy time mean [x_patch] of Lpatches - mean [x_patch] of Spatches
  ]

  ; Save timelapse (snapshot) image / Video frame
  if save_timelapse_img? and time >= next_tlapse_time [ save_tlapse_img ]
  if record_vid?   [ if ticks mod vid_rec_intval = 0 [carefully[ vid:record-view ][ print "Video capture failed."] ] ]

  ; All the main functions are below:
  unbind
  bind
  convert
  move

  tick
  set time (time + timestep);
end


to save_tlapse_img
  let inttime precision time 0
  let nZeros 3 - (length (word inttime))
  let $3digit_time ""
  repeat nZeros [set $3digit_time insert-item 0 $3digit_time "0"]
  set $3digit_time (word $3digit_time inttime)
  set-current-directory save-dir-name
  export-view (word run-index " t" $3digit_time ".png")
;  ifelse simple-savename? [    export-view (word run-index " t" $3digit_time ".png")      ]
;                          [    export-view (word file-prefix " t" $3digit_time ".png")    ]
  set next_tlapse_time    next_tlapse_time + tlapse_interval
end


to setup_world_from_input_file
  ; It's logically and computationally better to use "ifelse  - and what follows" rather than a series of "if's". However, doing so in Netlogo, the readability is greatly sacrificed.
  ; Therefore when the speed is not compromised much, I will use a series of "if's" instead of "ifelse - and what follows".
  if geometry-setup = "None" [
    set inpatches patches
    set outpatches no-patches
    ask inpatches [ set pcolor grey ]
  ]
  if geometry-setup = "Confinement" [
    import-pcolors-rgb input-geometry-fname
    set inpatches patches with [pcolor != [0 0 0]] ; when loading from image files, syntax like "black" does not work. Have to use RGB.
    set outpatches patches with [pcolor = [0 0 0]] ; when loading from image files, syntax like "black" does not work. Have to use RGB.
  ]
end


to add-sinusoidal-pip-perturbation
  if perturb-amplitude = 0 [user-message "Change perturb-amplitude to non-zero value"]
  if pert-wavelength = 0 [user-message "Change pert-wavelength (perturbation wavelength) ato non-zero value"]
  carefully [
    ask inpatches [
      let patch-xcor (pxcor / nGrid * worldLength)
      ;    let patch-ycor (pycor / nGrid * worldLength)
      let perturb_x (perturb-amplitude * cos (2 * 180 / pert-wavelength * patch-xcor))
      set x_patch (x_patch + perturb_x)
      if x_patch > 1 [ set x_patch 1 ]
      if x_patch < 0 [ set x_patch 0 ]
    ]
    ask inpatches [ represent-x-as-patch-color ]
    display
  ]
  [   user-message "Error. Change pert-wavelength (perturbation wavelength) to a non-zero value"   ]
end


to fluctuate
  let shot-noise-multiplier 1
  ask inpatches [
    ; Fluctuate pips
    let N_tot_pip  55555 * patchLength ^ 2
    let delta_xpatch random-normal 0      x_patch * (1 - x_patch) / sqrt(N_tot_pip)
    set x_patch x_patch + delta_xpatch * shot-noise-multiplier
    if x_patch < 0 [ set x_patch 0 ]
    if x_patch > 1 [ set x_patch 1 ]

    ; Fluctuate enzymes
    let N_kinase   N_tot_pip * k_mkon / k_koff
    let p_k patch_k_density / (k_mkon / k_koff)
    let q_k    1 - p_k
    let delta_kpatch random-normal 0   p_k * q_k / sqrt(N_kinase)
    set patch_k_density patch_k_density + delta_kpatch
    if patch_k_density < 0 [ set patch_k_density 0 ]

    let N_pptase    N_tot_pip * p_mkon / p_koff
    let p_p patch_p_density / (p_mkon / p_koff)
    let q_p    1 - p_p
    let delta_ppatch random-normal 0    p_p * q_p / sqrt(N_pptase)
    set patch_p_density patch_p_density + delta_ppatch
    if patch_p_density < 0 [ set patch_p_density 0 ]
  ]
  ask inpatches [ represent-x-as-patch-color ]
  display
end


to bind
  if Calculation-Type = "deterministic" [
    ask inpatches
    [
      ; kinase binding
      let k_binding-flux (k_mkon * (x_patch) * timestep)
      set patch_k_density (patch_k_density + k_binding-flux)
      ; phohsphatase binding
      if Enzyme-Pair-Type = "memK-memP"  [
        let p_binding-flux (p_mkon * (1 - x_patch) * timestep)
        set patch_p_density (patch_p_density + p_binding-flux)
      ]
    ]
  ]

  if Calculation-Type = "stochastic" [
    ask inpatches    [
      set k_Pon (k_mkon * x_patch * patchLength ^ 2 * timestep)
      if k_Pon > 1 [ user-message (word "k_Pon: " k_Pon) ]
      if random-float 1 < k_Pon
      [ sprout-kinases 1
        [
          set shape "circle"
          set size enz_size
          set color RGB-kinase
          if not show_enz?
          [ hide-turtle ]
      ]  ]

      if Enzyme-Pair-Type = "memK-memP"
      [
        set p_Pon (p_mkon * (1 - x_patch) * patchLength ^ 2 * timestep)
        if p_Pon > 1 [    user-message (word "p_Pon: " p_Pon)  ]
        if random-float 1 < p_Pon
        [ sprout-pptases 1
          [
            set shape "circle"
            set size enz_size
            set color RGB-pptase
            if not show_enz?
            [ hide-turtle ]
        ]  ]
      ]
    ]
  ]
end


to move
  ; Diffuse pip
  ask inpatches  [ set   updated_xpatch         x_patch * (1 - (4 * alpha-pip * (n_neighbors / 4))) + (sum [ x_patch ] of real_neighbors) * alpha-pip ]
  ask inpatches  [ set   x_patch                updated_xpatch  ]

  if Calculation-Type = "deterministic" [
    ask inpatches  [ set   updated_kpatch         patch_k_density * (1 - (4 * alpha-enz * (n_neighbors / 4))) + (sum [ patch_k_density ] of real_neighbors) * alpha-enz ]
    ask inpatches  [ set   patch_k_density                updated_kpatch  ]

    if Enzyme-Pair-Type = "memK-memP" [
      ask inpatches  [ set   updated_ppatch         patch_p_density * (1 - (4 * alpha-enz * (n_neighbors / 4))) + (sum [ patch_p_density ] of real_neighbors) * alpha-enz ]
      ask inpatches  [ set   patch_p_density                updated_ppatch  ]
    ]
  ]
  ; Move enzyems
  if Calculation-Type = "stochastic" [
    ask kinases
    [ if ( n_neighbors > 0 )    [     if random-float 1 < 1 - Pstay [ move-to one-of real_neighbors]     ]  ]
    ask pptases
    [ if ( n_neighbors > 0 )    [     if random-float 1 < 1 - Pstay [ move-to one-of real_neighbors]     ]  ]
  ]
end


to convert
  ask inpatches [

    let Kinase_contribution 0
    let Pptase_contribution 0
    let kinase-density 0
    let pptase-density 0

    if Calculation-Type = "deterministic" [
      set kinase-density patch_k_density
      set pptase-density patch_p_density
    ]

    if Calculation-Type = "stochastic" [
      set kinase-density count kinases-here / patchLength ^ 2
      set pptase-density count pptases-here / patchLength ^ 2
    ]

    ; Set the kinase's catalytic rate
    set Kinase_contribution   k_mkcat * kinase-density * (1 - x_patch) / (k_mKm + (1 - x_patch))

    ; Set the pptase's catalytic rate - The below code didn't use the "ifelse" syntax for better readability.
    if Enzyme-Pair-Type = "memK-memP" [ set Pptase_contribution   -1 * memP_mkcat * pptase-density * x_patch / (p_mKm + x_patch) ]
    if Enzyme-Pair-Type = "memK-solP" [ set Pptase_contribution   -1 * solP_mkcat * x_patch / (p_mKm + x_patch) ]

    set dx_patch (Kinase_contribution + Pptase_contribution) * timestep

    if disallow-too-large-dx? and abs(dx_patch) > 0.5 [ user-message (word "abs(d_x_patch) > 0.5 patch: " dx_patch " at " pxcor " " pycor)  ]
    set x_patch (x_patch + dx_patch)
    if x_patch < 0 [ set x_patch 0 ]
    if x_patch > 1 [ set x_patch 1 ]
  ]
end


to unbind
  if Calculation-Type = "deterministic" [
    ask inpatches [
      let k_unbinding-flux k_koff * patch_k_density * timestep
      set patch_k_density patch_k_density - k_unbinding-flux
      if patch_k_density < 0 [ set patch_k_density 0 ]

      if Enzyme-Pair-Type = "memK-memP" [
        let p_unbinding-flux p_koff * patch_p_density * timestep
        set patch_p_density patch_p_density - p_unbinding-flux
        if patch_p_density < 0 [ set patch_p_density 0 ]
      ]
    ]
  ]

  if Calculation-Type = "stochastic" [
    ask kinases [ if random-float 1 < const-k_Poff [ die ]  ]
    if Enzyme-Pair-Type = "memK-memP" [   ask pptases  [ if random-float 1 < const-p_Poff [ die ]  ]    ]
  ]
end


to set-neighbors_and_outColors
  ask inpatches
  [ set real_neighbors neighbors4 with [member? self inpatches]
    set n_neighbors count real_neighbors  ]
  ; Color the outpatches brown, in order to make it (visually) obvious that this is simulation ;
  if geometry-setup = "Confinement"
  [ ask outpatches
    [ set pcolor colorname-outpatches]
    ask outpatches
    [ if any? neighbors with [member? self inpatches] [set pcolor black] ] ]
end


to represent-x-as-patch-color
  let x1 map [[x] -> x_patch * x ] RGB-pip2
  let x2 map [[x] -> (1 - x_patch) * x ] RGB-pip1 ;  If 4 components vs. 3 components matter is involved, consider this code:  let x2 map [[x] -> (1 - x_initial) * x ] but-last RGB-pip1
  set pcolor (map + x1 x2)
end


to-report retrieve-simul_info_string
  if Enzyme-Pair-Type = "memK-memP" [
    report (word (substring Calculation-Type 0 3) " " worldLength "um "
      "mkon-off-mkcat-mKm of memK - " k_mkon "-" k_koff "-" k_mkcat "-" k_mKm
      "_ of memP - " p_mkon "-" p_koff "-" memP_mkcat "-" p_mKm    )  ]

  if Enzyme-Pair-Type = "memK-solP" [
    report (word (substring Calculation-Type 0 3) " " worldLength "um "
       "on-off-cat-Km of memK - " k_mkon "-" k_koff "-" k_mkcat "-" k_mKm
      "_ of solP - " solP_mkcat "-" p_mKm    )  ]
end


;to-report retrieve-num-date-time  ; current date in numerical format, yyyy-mm-dd
;  let $dt substring date-and-time 16 27  ; "21-May-2013"
;  let $dt2 substring date-and-time 0 8  ; 01:19
;  report (word (substring $dt 7 11)           ; yyyy
;    "-" (retrieve-month-num substring $dt 3 6)  ; mm
;    "-" (substring $dt 0 2)           ; dd
;    "_" (substring $dt2 0 2) "_" (substring $dt2 3 5) "_" (substring $dt2 6 8))
;end
;
;
;to-report retrieve-month-num [ #mon ]
;  let $index 1 + position #mon
;    ["Jan""Feb""Mar""Apr""May""Jun""Jul""Aug""Sep""Oct""Nov""Dec"]
;  report substring (word (100 + $index)) 1 3  ; force 2-digit string
;end


; Modules not being used at the moment
to click-x-up
  if mouse-down? and (member? patch mouse-xcor mouse-ycor inpatches)  [
    ask patch mouse-xcor mouse-ycor [
      ask patches in-radius 3 [set x_patch 1 ]
    ]
  ]
  ask inpatches [ represent-x-as-patch-color ]
  display
;  tick
end
@#$#@#$#@
GRAPHICS-WINDOW
11
375
419
784
-1
-1
8.0
1
15
1
1
1
0
1
1
1
0
49
0
49
1
1
1
ticks
30.0

BUTTON
16
35
125
75
NIL
setup
NIL
1
T
OBSERVER
NIL
R
NIL
NIL
1

BUTTON
15
131
126
171
go
if run-index < N-runs [go]
T
1
T
OBSERVER
NIL
G
NIL
NIL
1

PLOT
1049
360
1249
573
PIP fraction
time
NIL
0.0
0.0
0.0
0.0
true
true
"" ""
PENS
"pip1" 1.0 2 -13345367 true "" "plotxy time 1 - avg_x"
"pip2" 1.0 2 -4079321 true "" "plotxy time avg_x"

PLOT
1052
91
1251
305
Number of Enzymes
time
NIL
0.0
0.0
0.0
0.0
true
true
"" ""
PENS
"K" 1.0 2 -955883 true "" "if Calculation-Type = \"deterministic\" [plotxy time sum ([patch_k_density] of inpatches) * patchLength ^ 2]\nif Calculation-Type = \"stochastic\" [plotxy time count kinases]\n"
"P" 1.0 2 -11221820 true "" "if Calculation-Type = \"deterministic\" [plotxy time sum ([patch_p_density] of inpatches) * patchLength ^ 2]\nif Calculation-Type = \"stochastic\" [plotxy time count pptases]\n"

INPUTBOX
530
414
618
474
k_mkon
0.1
1
0
Number

INPUTBOX
626
414
715
474
k_koff
0.7
1
0
Number

INPUTBOX
723
414
813
474
k_mkcat
10.0
1
0
Number

INPUTBOX
531
528
619
588
p_mkon
0.02
1
0
Number

INPUTBOX
625
527
714
587
p_koff
0.1
1
0
Number

INPUTBOX
723
491
815
551
memP_mkcat
15.0
1
0
Number

SWITCH
24
331
161
364
show_enz?
show_enz?
0
1
-1000

INPUTBOX
627
320
714
380
D_pip
2.0
1
0
Number

SLIDER
172
331
415
364
enz_size
enz_size
0
1
0.49
0.01
1
NIL
HORIZONTAL

SLIDER
269
286
414
319
world_pixel_length
world_pixel_length
100
700
400.0
50
1
NIL
HORIZONTAL

PLOT
1051
581
1250
795
Max Pon-patch
time
NIL
0.0
0.0
0.0
0.0
true
true
"" ""
PENS
"K" 1.0 0 -3844592 true "" "plotxy time max [k_Pon] of patches"
"P" 1.0 0 -13403783 true "" "plotxy time max [p_Pon] of patches"

CHOOSER
748
334
840
379
nGrid
nGrid
1 2 3 4 5 6 7 9 10 12 15 16 18 20 24 25 27 30 34 36 40 41 50 51 59 60 64 66 70 89 90 96 100 150 200 250 350 500
22

SWITCH
534
664
706
697
save_timelapse_img?
save_timelapse_img?
0
1
-1000

SWITCH
733
691
853
724
record_vid?
record_vid?
0
1
-1000

SWITCH
1053
52
1469
85
save_all_plots?
save_all_plots?
1
1
-1000

INPUTBOX
535
224
626
284
endtime
2.0
1
0
Number

INPUTBOX
376
32
805
102
input-geometry-fname
C:\\Users\\Neil\\Dropbox\\Research\\20200824\\confinements\\500-snail6.png
1
0
String

INPUTBOX
845
319
923
379
worldLength
30.0
1
0
Number

INPUTBOX
630
224
714
284
timestep
0.01
1
0
Number

CHOOSER
195
35
366
80
geometry-setup
geometry-setup
"None" "Confinement"
1

SLIDER
731
730
853
763
vid_rec_intval
vid_rec_intval
10
300
10.0
10
1
NIL
HORIZONTAL

SLIDER
534
703
706
736
tlapse_interval
tlapse_interval
0.5
150
0.5
0.5
1
s
HORIZONTAL

SWITCH
1269
281
1467
314
disallow-too-large-dx?
disallow-too-large-dx?
0
1
-1000

PLOT
1267
322
1467
533
Max Min dx_patch
time
NIL
0.0
0.0
0.0
0.0
true
false
"" ""
PENS
"max" 1.0 0 -16777216 true "" "if inpatches != 0 [plotxy time max [dx_patch] of inpatches ]"
"min" 1.0 0 -7500403 true "" "if inpatches != 0 [plotxy time min [dx_patch] of inpatches ]"

MONITOR
927
334
1003
379
patchLength
worldLength / nGrid
7
1
11

SWITCH
24
287
197
320
timestamp-on-image?
timestamp-on-image?
0
1
-1000

INPUTBOX
820
414
910
474
k_mKm
2.0
1
0
Number

INPUTBOX
821
527
908
587
p_mKm
0.5
1
0
Number

INPUTBOX
723
556
816
616
solP_mkcat
15.0
1
0
Number

INPUTBOX
531
320
620
380
D_enz
0.2
1
0
Number

INPUTBOX
1606
99
1706
159
pert-wavelength
1.0
1
0
Number

BUTTON
1503
59
1705
92
add-sinusoidal-pip-perturbation
add-sinusoidal-pip-perturbation
NIL
1
T
OBSERVER
NIL
Q
NIL
NIL
1

INPUTBOX
1502
98
1599
158
perturb-amplitude
1.0
1
0
Number

BUTTON
196
87
365
120
change input geometry file
carefully [ set input-geometry-fname user-file ]\n  [ print \"geometry-input-fname not updated.\" ]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
21
209
213
254
Calculation-Type
Calculation-Type
"stochastic" "deterministic"
0

CHOOSER
223
208
415
253
Enzyme-Pair-Type
Enzyme-Pair-Type
"memK-memP" "memK-solP"
1

INPUTBOX
847
63
1002
123
KIN-PPT-X
[.0 .0 .5]
1
0
String

BUTTON
847
128
1001
161
default value
set KIN-PPT-X \"[.0 .0 .5]\"
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
1503
245
1704
278
smaller-part-test?
smaller-part-test?
1
1
-1000

MONITOR
1502
162
1704
207
pert-angular wavenumber
2 * pi / pert-wavelength
5
1
11

MONITOR
205
281
262
326
NIL
time
3
1
11

PLOT
1267
91
1466
267
dxdt vs x
NIL
NIL
0.0
1.0
0.0
0.0
true
false
"" ""
PENS
"default" 1.0 0 -14835848 true "" ""
"dxdt0" 1.0 0 -16777216 true "" ""
"x0.5" 1.0 0 -16777216 true "" ""

CHOOSER
750
231
843
276
N-runs
N-runs
1 2 3 4 5 9 10 20 30 40 50
2

MONITOR
850
232
928
277
Current run
run-index + 1
17
1
11

PLOT
1266
582
1466
795
xL-xS
time
NIL
0.0
0.1
-1.1
1.1
true
false
"" ""
PENS
"0" 1.0 2 -7500403 true "" ""
"1" 1.0 2 -2674135 true "" ""
"2" 1.0 2 -13403783 true "" ""
"3" 1.0 2 -7171555 true "" ""
"4" 1.0 0 -14070903 true "" ""
"5" 1.0 2 -2064490 true "" ""
"6" 1.0 2 -13840069 true "" ""
"7" 1.0 2 -9276814 true "" ""
"8" 1.0 2 -11221820 true "" ""
"9" 1.0 2 -16777216 true "" ""

BUTTON
196
135
366
168
set-save-dir
carefully [\nset save-dir-name user-directory\nset setup-success? false\n]\n[ \nuser-message (word \"Save directory unchanged.\")\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
1266
545
1465
578
plot-xL-xS?
plot-xL-xS?
0
1
-1000

SWITCH
534
744
706
777
simple-savename?
simple-savename?
1
1
-1000

MONITOR
1058
310
1150
355
avg  k-density
count kinases / (patchLength ^ 2 * count inpatches)
5
1
11

MONITOR
1155
310
1248
355
avg  p-density
count pptases / (patchLength ^ 2 * count inpatches)
5
1
11

TEXTBOX
1508
34
1714
69
Deterministic Simulation Settings
13
0.0
1

TEXTBOX
537
297
687
315
Diffusion constants
13
0.0
1

TEXTBOX
534
391
684
409
Kinase parameters
13
0.0
1

TEXTBOX
533
492
683
510
Phosphatase parameters
13
0.0
1

TEXTBOX
1056
32
1206
50
Plots
13
0.0
1

TEXTBOX
536
637
729
656
Saving Images and Videos
13
0.0
1

TEXTBOX
754
294
958
326
Membrane size and patches
13
0.0
1

INPUTBOX
377
104
806
168
save-dir-name
C:\\Users\\Neil\\Dropbox\\github\\PIPpolarize\\results\\123\\
1
0
String

MONITOR
15
81
125
126
NIL
setup-success?
17
1
11

TEXTBOX
539
203
689
221
Time setting
13
0.0
1

TEXTBOX
23
186
173
204
Simulation Type Setting
13
0.0
1

TEXTBOX
753
208
903
226
Simulation Runs
13
0.0
1

TEXTBOX
24
267
174
285
Visual Settings
13
0.0
1

TEXTBOX
850
37
1000
55
Starting condition
13
0.0
1

TEXTBOX
1506
223
1656
241
Special Settings\t
13
0.0
1

TEXTBOX
21
10
171
28
Main Commands\t
13
0.0
1

TEXTBOX
201
11
351
29
File I/O
13
0.0
1

SWITCH
1266
800
1466
833
save-xL-xS?
save-xL-xS?
0
1
-1000

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="COMB" repetitions="10" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>reset</setup>
    <go>go</go>
    <final>set simple-save-name? false</final>
    <enumeratedValueSet variable="setsavedir?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simple-save-name?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Enzyme-Pair-Type">
      <value value="&quot;mem K - mem P&quot;"/>
      <value value="&quot;mem K - sol P&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="COMB" repetitions="9" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>reset</setup>
    <go>go</go>
    <final>set simple-save-name? false</final>
    <enumeratedValueSet variable="setsavedir?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simple-save-name?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
