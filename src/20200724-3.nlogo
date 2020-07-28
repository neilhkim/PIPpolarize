extensions [ vid ]

globals [
  ; Internally set Constants ;
  nlogoColor-outpatches
  RGB-pip1
  RGB-pip2
  RGB-kinase
  RGB-pptase

  ; Calculated Constants ;
  timestep
  dist_enz ; distance travelled by enzyme per timestep
  pip-diffuse-out-portion ; how much PIP get dispersed to neighbors
  enz-diffuse-out-portion ; (when assuming continuum) how much enzymes get dispersed to neighbors
  patchLength
  k_Poff
  p_Poff

  ; Patch lists ;
  inpatches
  outpatches
  edgeInPatches
  non_Edge_Patches

  ; Trackers - ever changing ;
  avg_x
  time
  next_snapshot_time

  file-prefix
]


breed [ pip1s pip1 ]
breed [ pip2s pip2 ]
breed [ kinases kinase ]
breed [ pptases pptase ]
breed [ clocks clock ]


turtles-own [
  orig_heading lr_collision? ud_collision? orig_dx orig_dy
  mydist
  ;k_or_p - If you don't want kinase always overlapping on top of pptase, then try: assigning them to the same breed, but internally distinguishing between K and P.
]


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


to reset
  clear-all
  ; Fix constants (internally defined) ;
  set nlogoColor-outpatches brown
  set RGB-pip1 [0 100 255] ;  set RGB-pip1 extract-rgb blue
  set RGB-pip2 [255 200 0]
  set RGB-kinase [255 150 0]
;  set RGB-pptase [0 0 200]
  set RGB-pptase [0 150 255]

  ; Setup space ;
  set inpatches no-patches
  set outpatches no-patches
  set edgeInPatches no-patches
  set non_Edge_Patches no-patches

  ; Space - Reflect user-input
  resize-world 0 (nGrid - 1) 0 (nGrid - 1)
  if ten-percent-height-test? [resize-world 0 (nGrid - 1) 0 (nGrid - 1) / 10]
  set-patch-size world_pixel_length / nGrid
  ifelse wrap? [ __change-topology true  true  ]
                [ __change-topology false false ]
  setup_world_from_input_file
  set patchLength worldLength / nGrid
  set_neighbors_and_outColors
  setup_edgepatches

  RESET-TICKS
  if d_t-input-option = "exponent"
  [ set timestep 10 ^ log_timestep
    set timestep_input 0
  ]
  if d_t-input-option = "direct"
  [ set timestep timestep_input
    set log_timestep 999   ]

  ; Set-up clock
  if display_time?
  [  create-clocks 1 [ set label (precision time 1) setxy (min-pxcor + nGrid / 5) (max-pycor - nGrid / 20) set size 0 ]  ]
  ;crt 1 [ set label (word "q=" (precision (2 * pi / lambda) 4)) setxy (max-pxcor - nGrid / 5) (max-pycor - nGrid / 20) set size 0 ]

  ; Setup (unchanging) rates
  set pip-diffuse-out-portion 1 - exp(-4 * D_pip * timestep / (worldLength / nGrid) ^ 2)
  if pip-diffuse-out-portion > 1 [ user-message "PIP diffusion fraction is more than 1" ]
  set enz-diffuse-out-portion 1 - exp(-4 * D_enz * timestep / (worldLength / nGrid) ^ 2)
  if enz-diffuse-out-portion > 1 [ user-message "Enzyme diffusion fraction is more than 1" ]
  set dist_enz sqrt(4 * D_enz * timestep) * (1 / worldLength * nGrid)
  set k_Poff k_koff * timestep
  set p_Poff p_koff * timestep
  if k_Poff > 1 [user-message "k_off_prob_sto > 1"]
  if p_Poff > 1 [user-message "p_off_prob_sto > 1"]


  set file-prefix retrieve-simul_info_string
  if setsavedir? [initialize-file-saving]
  if record_vid?  [ carefully [vid:start-recorder][ print error-message ] ]
  initialize_patches
  ; Initialize times
  set time 0

  plot_dxdt_vs_x

end

to plot_dxdt_vs_x
  let Nxpoints 100
  let xpoints n-values (Nxpoints + 1) [ x -> x / 100 ]
  set-current-plot "dxdt vs x"



  if Enzyme-Pair-Type = "mem K - sol P" [
    foreach xpoints [
      x -> let dxdt 0
      if k_S-dependence = "linear" [ set dxdt k_kon / k_koff * x * k_kcat * (1 - x) ]
      if k_S-dependence = "MM" [ set dxdt k_kon / k_koff * x * k_kcat * (1 - x) / (k_Km + 1 - x) ]
      if p_S-dependence = "linear" [ set dxdt dxdt - psol_kcat_per_A * x ]
      if p_S-dependence = "MM" [ set dxdt dxdt - psol_kcat_per_A * x / (p_Km + x) ]
      plotxy x dxdt
    ]
    set-current-plot-pen "dxdt0"
    foreach xpoints [x -> plotxy x 0]
  ]
  if Enzyme-Pair-Type = "mem K - mem P" [
    foreach xpoints [
      x -> let dxdt 0
      if k_S-dependence = "linear" [ set dxdt k_kon / k_koff * x * k_kcat * (1 - x) ]
      if k_S-dependence = "MM" [ set dxdt k_kon / k_koff * x * k_kcat * (1 - x) / (k_Km + 1 - x) ]
      if p_S-dependence = "linear" [ set dxdt dxdt - p_kon / p_koff * (1 - x) * p_kcat * x ]
      if p_S-dependence = "MM" [ set dxdt dxdt - p_kon / p_koff * (1 - x) * p_kcat * x / (p_Km + x) ]
      plotxy x dxdt
    ]
    set-current-plot-pen "dxdt0"
    foreach xpoints [x -> plotxy x 0]
  ]
  set-current-plot-pen "x0.5"
  plotxy 0.5 plot-y-max
  plotxy 0.5 plot-y-min
end

to initialize_patches
    ; Initialized Patches
  let init-k_patch-density item 0 read-from-string string-KPX
  let init-p_patch-density item 1 read-from-string string-KPX
  if Enzyme-Pair-Type = "mem K - sol P" and init-p_patch-density != 0
  [ print "**Warning: \"solution Pptase setup\" and still adding binding Pptase on the membrane as initial condition?"]
  let init-x_patch item 2 read-from-string string-KPX

  ask kinases [die]
  ask pptases [die]
  ask inpatches [
    set patch_k_density 0
    set patch_p_density 0
  ]

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
        [ set mydist dist_enz
          set lr_collision? false
          set ud_collision? false
          let exit_flag false
          while [ exit_flag = false ]
          [ setxy random-xcor random-ycor
            if member? patch-here inpatches [ set exit_flag true ]
          ]
          set shape "circle"
          set size enz_size
          set color RGB-kinase
          if not show_enz?
          [ hide-turtle ]
    ]  ]  ]
    let n_init-pptases init-p_patch-density * count inpatches * patchLength ^ 2
    repeat n_init-pptases[
      ask one-of inpatches [
        sprout-pptases 1
        [ set mydist dist_enz
          set lr_collision? false
          set ud_collision? false
          let exit_flag false
          while [ exit_flag = false ]
          [ setxy random-xcor random-ycor
            if member? patch-here inpatches [ set exit_flag true ]
          ]
          set shape "circle"
          set size enz_size
          set color RGB-kinase
          if not show_enz?
          [ hide-turtle ]
    ]  ]  ]
  ]
  ask inpatches [   set x_patch init-x_patch   ]
  if initial-fluctuation? or Calculation-Type = "deterministic" [ fluctuate ]
  ask inpatches [ update_pcolor ]
  set avg_x mean [x_patch] of inpatches
end


to go
  if display_time? [ ask clocks [ set label (precision time 2) ] ]

  ; Snapshot condition check
  if save_timelapse_img? and time >= next_snapshot_time [
    let inttime precision time 0
    let nfill 3 - (length (word inttime))
    let $3digit_time ""
    repeat nfill [set $3digit_time insert-item 0 $3digit_time "0"]
    set $3digit_time (word $3digit_time inttime)

    if snap-w-and-wo-enz?
    [
      ask kinases [show-turtle]
      ask pptases [show-turtle]
      set-patch-size patch-size * 10
      display
      ;export-view (word "ENZ" retrieve-simul_info_string " t" $3digit_time ".png")
      export-view (word "ENZ" file-prefix " t" $3digit_time ".png")
      set-patch-size patch-size / 10
      ask kinases [hide-turtle]
      ask pptases [hide-turtle]
      display
    ]
    ;export-view (word retrieve-simul_info_string " t" $3digit_time ".png")

    ifelse simple-save-name? [
      let datetime substring file-prefix 0 19
      export-view (word datetime " " $3digit_time ".png")
    ][
      export-view (word file-prefix " t" $3digit_time ".png")
    ]
    set next_snapshot_time next_snapshot_time + tlapse_interval
  ]

  if record_vid?
  [ if ticks mod vid_rec_intval = 0 [ vid:record-view ] ]

  if show_enz? [ ask kinases [show-turtle] ask pptases [show-turtle] ]
  if not show_enz? [ ask kinases [hide-turtle] ask pptases [hide-turtle] ]

  ; All the main functions are below:
  unbind
  convert
  move
  bind
  ask inpatches [ update_pcolor ]
  set avg_x mean [x_patch] of inpatches

  if (avg_x < loss_threshold) or (avg_x > 1 - loss_threshold) or (time > endtime)
  [
    ;if record_vid? [  vid:save-recording (word retrieve-simul_info_string "_mov.mp4") ]
    if record_vid? [  vid:save-recording (word file-prefix "_mov.mp4") ]
    ;if save_timelapse_img? or record_vid? or save_plots?  [export-interface (word "iface-End " retrieve-simul_info_string ".png")]
    if save_timelapse_img? or record_vid? or save_plots?  [export-interface (word "iface-End " file-prefix ".png")]
    if save_plots? [export-all-plots "allplots.csv"]
    stop
  ]

; Clear non-accumulative plots
;  set-current-plot "Max Pon-patch" clear-plot
;  set-current-plot "mean_x of ROIs" clear-plot
;  set-current-plot "Kinase count - Edge vs non-Edge" clear-plot
  ; Finish-up Recording

  tick
  set time (time + timestep);
end


to change-input-geometry-file
  carefully [ set input-geometry-fname user-file ]
  [ print "geometry-input-fname not updated." ]
end

to setup_world_from_input_file
  ; It's logically and computationally better to use "ifelse  - and what follows" rather than a series of "if's". However, doing so in Netlogo, the readability is greatly sacrificed.
  ; Therefore when the speed is not compromised very much, I will use a series of "if's" instead of "ifelse - and what follows".
  if geometry-setup = "None" [
    set inpatches patches
    set outpatches no-patches
    ask inpatches [ set pcolor grey ]
  ]
  if geometry-setup = "Confinement" [
    import-pcolors-rgb input-geometry-fname
    set inpatches patches with [pcolor != [0 0 0]] ; loading from image files, syntax like "black" does not work. Have to use RGB.
    set outpatches patches with [pcolor = [0 0 0]] ; loading from image files, syntax like "black" does not work. Have to use RGB.
    ; Checking if enzymes are sticking or excluded from the boundary patches
  ]
end

to setup_edgepatches
  set edgeInPatches inpatches with [n_neighbors < 8];
  set non_Edge_Patches inpatches with [not member? self edgeInPatches]
end


to perturb-pips
  ask inpatches [

    let x (pxcor / nGrid * worldLength)
;    print x
    let y (pycor / nGrid * worldLength)
;    let perturb_x (perturb-amplitude * cos (180 / 2 / pi * q * x) * cos (180 / 2 / pi * q * y) )
    let perturb_x (perturb-amplitude * cos (2 * 180 / lambda * x))
;    let perturb_x (perturb-amplitude * cos (2 * pi * q * x))
    set x_patch (x_patch + perturb_x)



    if x_patch > 1 [ set x_patch 1 ]
    if x_patch < 0 [ set x_patch 0 ]
  ]
  ask inpatches [ update_pcolor ]
  display
end


to fluctuate
  ask inpatches [
    let N 55555
    let delta_xpatch random-normal 0 (x_patch * (1 - x_patch) / sqrt(N * patchLength ^ 2))
;    set delta_xpatch delta_xpatch * 5000

;    if delta_xpatch > 0 [ print delta_xpatch ]
    set x_patch x_patch + delta_xpatch ;* 1000
    if x_patch < 0 [ set x_patch 0 ]
    if x_patch > 1 [ set x_patch 1 ]

    let KN N / 100
    let p_k patch_k_density / (k_kon / k_koff)
    let q_k 1 - p_k
    let delta_kpatch random-normal 0 (p_k * q_k / sqrt(KN * patchLength ^ 2))
    set patch_k_density patch_k_density + delta_kpatch
    if patch_k_density < 0 [ set patch_k_density 0 ]

    let KP N / 100
    let p_p patch_p_density / (p_kon / p_koff)
    let q_p 1 - p_p
    let delta_ppatch random-normal 0 (p_p * q_p / sqrt(KP * patchLength ^ 2))
    set patch_p_density patch_p_density + delta_ppatch
    if patch_p_density < 0 [ set patch_p_density 0 ]
  ]
  ask inpatches [ update_pcolor ]
  display
end

;to fluctuate_pips
;  ask inpatches [
;    let N 55555
;    let delta_xpatch random-normal 0 (x_patch * (1 - x_patch) / sqrt(N * patchLength ^ 2))
;;    if delta_xpatch > 0 [ print delta_xpatch ]
;    set x_patch x_patch + delta_xpatch
;    if x_patch < 0 [ set x_patch 0 ]
;    if x_patch > 1 [ set x_patch 1 ]
;  ]
;  ask inpatches [ update_pcolor ]
;  display
;end


to set-KPT-default
  set string-KPX "[.0 .0 .5]"
end

to bind
  if Calculation-Type = "deterministic" [
    let k_binding-flux 0
    let p_binding-flux 0
    ask inpatches
    [
      set k_binding-flux (k_kon * (x_patch) * timestep)
      set patch_k_density (patch_k_density + k_binding-flux)

      if Enzyme-Pair-Type = "mem K - mem P"  [
        set p_binding-flux (p_kon * (1 - x_patch) * timestep)
        set patch_p_density (patch_p_density + p_binding-flux)
      ]
    ]
  ]

  if Calculation-Type = "stochastic" [
    ask inpatches
    [
      set k_Pon (k_kon * x_patch * patchLength ^ 2 * timestep)
      if k_Pon > 1 [ user-message (word "k_Pon: " k_Pon) ]
      if random-float 1 < k_Pon
      [ sprout-kinases 1
        [ set mydist dist_enz
          set lr_collision? false
          set ud_collision? false
          setxy xcor + random-float 1 - 0.5 ycor + random-float 1 - 0.5
          set shape "circle"
          set size enz_size
          set color RGB-kinase
          if not show_enz?
          [ hide-turtle ]
      ]  ]

      if Enzyme-Pair-Type = "mem K - mem P"
      [
        set p_Pon (p_kon * (1 - x_patch) * patchLength ^ 2 * timestep)
        if p_Pon > 1 [    user-message (word "p_Pon: " p_Pon)  ]
        if random-float 1 < p_Pon
        [ sprout-pptases 1
          [ set mydist dist_enz
            set lr_collision? false
            set ud_collision? false
            setxy xcor + random-float 1 - 0.5 ycor + random-float 1 - 0.5
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
  ask inpatches  [ set   updated_xpatch         x_patch * (1 - pip-diffuse-out-portion * (n_neighbors / 8)) + (sum [ x_patch ] of real_neighbors) / 8 * pip-diffuse-out-portion  ]
  ask inpatches  [ set   x_patch                updated_xpatch  ]

  if Calculation-Type = "deterministic" [
    ask inpatches  [  set updated_kpatch patch_k_density * (1 - enz-diffuse-out-portion * (n_neighbors / 8)) + (sum [ patch_k_density ] of real_neighbors) / 8 * enz-diffuse-out-portion  ]
    ask inpatches  [  set patch_k_density updated_kpatch  ]

    if Enzyme-Pair-Type = "mem K - mem P" [
      ask inpatches    [ set updated_ppatch patch_p_density * (1 - enz-diffuse-out-portion * (n_neighbors / 8)) + (sum [ patch_p_density ] of real_neighbors) / 8 * enz-diffuse-out-portion ]
      ask inpatches    [ set patch_p_density updated_ppatch ]
    ]
  ]

  if Calculation-Type = "stochastic" [
    ; Move Enzymes (turtles)
    if Enzymes-touch-the-edge? = true
    [ ; Move carelessly (simply approach the boundary and stop when you're close)
      ask turtles ; This asks to clocks, too. But clocks won't move, so this is okay.
      [ rt random 360
        ifelse (patch-ahead mydist = nobody or not member? (patch-ahead mydist) inpatches)
        [ while [patch-ahead (mydist / 10) != nobody and member? (patch-ahead (mydist / 10)) inpatches] [fd mydist / 10 ] ] ; if
        [ fd mydist ]  ; else
      ]
    ]

    if Enzymes-touch-the-edge? = false
    [ ; Move carefully
      ask turtles
      [ right random 360
        set lr_collision? false
        set ud_collision? false

        if (member? (patch-ahead mydist) outpatches) or (patch-ahead mydist = nobody)
        [ set orig_heading heading
          set orig_dx dx
          set orig_dy dy
          ; Test collision to left/right (for left, dx will be negative)
          set heading 90
          if (member? (patch-ahead (orig_dx * mydist)) outpatches) or (patch-ahead (orig_dx * mydist) = nobody )
          [ set lr_collision? true ]
          ; Test collision to up/down (for down, dx will be negative)
          set heading 0
          if (member? (patch-ahead (orig_dy * mydist)) outpatches) or (patch-ahead (orig_dy * mydist) = nobody)
          [ set ud_collision? true ]
          ; There is a special case where only when moving in the original heading, it will meet an obstacle. ; It's when running to a protruded corner. Then, depending on the incident angle, the reflection direction is determined.
          set heading orig_heading
          if lr_collision? = false and ud_collision? = false and (member? (patch-ahead mydist) outpatches)
          [ let patchx [pxcor] of patch-ahead mydist
            let patchy [pycor] of patch-ahead mydist
            ifelse (patchy - 0.5 - ycor) / (patchx - 0.5 - xcor) > dy / dx
            [ set ud_collision? true ] ; if
            [ set lr_collision? true ] ; else
          ]
          ; Apply the above determined reflection
          set heading orig_heading
          if lr_collision?
          [ set heading (- heading)
            set lr_collision? false ]
          if ud_collision?
          [ set heading (180 - heading)
            set ud_collision? false ]  ]
        ; Now, you can finally go forward
        forward mydist

        ; Now there is a special case where doing all of the above still makes the enzyme go over the boundary. This is when 2+ reflections is needed. In this case, I say the timestep is too big and stop the run.
        if (member? patch-here outpatches)
        [ set heading heading + 180
          pen-down
          forward mydist
          user-message "Reduce timescale. Model detects more than 1 reflection per timestep."
          pen-erase ]
      ]
    ]
  ]
end


to convert

  ask inpatches [

    let Kinase_contribution 0
    let Pptase_contribution 0
    let general-kinase-density 0
    let general-pptase-density 0

    if Calculation-Type = "deterministic" [
      set general-kinase-density patch_k_density
      set general-pptase-density patch_p_density
    ]

    if Calculation-Type = "stochastic" [
      set general-kinase-density count kinases-here / patchLength ^ 2
      set general-pptase-density count pptases-here / patchLength ^ 2
    ]


    ; Set the kinase's catalytic rate
    ifelse k_S-dependence = "none"  [  set Kinase_contribution   k_kcat * general-kinase-density  ]
    [ ifelse k_S-dependence = "linear" [  set Kinase_contribution   k_kcat * general-kinase-density * (1 - x_patch)  ]
      [ ifelse k_S-dependence = "MM" [  set Kinase_contribution   k_kcat * general-kinase-density * (1 - x_patch) / (k_Km + (1 - x_patch))   ] [ print "K substrate-dependence setup error."  ]  ]
    ]

    ; Set the pptase's catalytic rate
    if Enzyme-Pair-Type = "mem K - mem P" [
      ifelse p_S-dependence = "none" [ set Pptase_contribution   -1 * p_kcat * general-pptase-density ]
      [ ifelse p_S-dependence = "linear" [ set Pptase_contribution   -1 * p_kcat * general-pptase-density * x_patch ]
        [ ifelse p_S-dependence = "MM" [ set Pptase_contribution   -1 * p_kcat * general-pptase-density * x_patch / (p_Km + x_patch) ] [ print "P substrate-dependence setup error." ]  ]
      ]
    ]
    if Enzyme-Pair-Type = "mem K - sol P" [
      ifelse p_S-dependence = "none" [ set Pptase_contribution   -1 * psol_kcat_per_A  ]
      [ ifelse p_S-dependence = "linear" [ set Pptase_contribution   -1 * psol_kcat_per_A * x_patch ]
        [ ifelse p_S-dependence = "MM" [ set Pptase_contribution   -1 * psol_kcat_per_A * x_patch / (p_Km + x_patch) ] [ print "P substrate-dependence setup error." ]  ]
      ]
    ]

    set dx_patch (Kinase_contribution + Pptase_contribution) * timestep

    if disallow-too-large-dx? and abs(dx_patch) > 0.5 [ user-message (word "abs(d_x_patch) > 0.5 patch: " dx_patch " at " pxcor " " pycor)  ]
    set x_patch (x_patch + dx_patch)
    if x_patch < 0 [ set x_patch 0 ]
    if x_patch > 1 [ set x_patch 1 ]

  ]
end


to unbind
  if Calculation-Type = "deterministic" [
    ask inpatches[
      let k_unbinding-flux k_koff * patch_k_density * timestep
      set patch_k_density patch_k_density - k_unbinding-flux
      if patch_k_density < 0 [ set patch_k_density 0 ]

      if Enzyme-Pair-Type = "mem K - mem P" [
        let p_unbinding-flux p_koff * patch_p_density * timestep
        set patch_p_density patch_p_density - p_unbinding-flux
        if patch_p_density < 0 [ set patch_p_density 0 ]
      ]
    ]
  ]

  if Calculation-Type = "stochastic" [
    ask kinases [ if random-float 1 < k_Poff [ die ]  ]
    if Enzyme-Pair-Type = "mem K - mem P" [
      ask pptases  [ if random-float 1 < p_Poff [ die ]  ]
    ]
  ]
end


to save-recorded
  ;vid:save-recording (word retrieve-simul_info_string "_conj_kin_vid.mp4")
  vid:save-recording (word file-prefix "_conj_kin_vid.mp4")
end


to initialize-file-saving
  if save_timelapse_img? or record_vid? or save_plots?  [
    carefully [
      set-current-directory user-directory
      export-interface (word "iface-t0 " file-prefix ".png")
    ]
    [ print "Results saving directory not selected." ] ]


end


to run_edge_test
  reset
  if count edgeInPatches = 0
  [ user-message "There are no edge patches. Check 'wrap' setting." ]

  create-kinases 1000
  [ set mydist dist_enz
    set lr_collision? false
    set ud_collision? false
    set shape "circle"
    set size enz_size
    set color RGB-kinase
    setxy random-xcor random-ycor
  ]

  while [time < endtime]
  [ tick
    set time (time + timestep)
    move ]
end


to set_neighbors_and_outColors
  ask inpatches
  [ set real_neighbors neighbors with [member? self inpatches]
    set n_neighbors count real_neighbors  ]
  ; Decorate the outpatches to make it visually obvious that this is simulation ;
  if geometry-setup = "Confinement"
  [ ask outpatches
    [ set pcolor nlogoColor-outpatches - 1 ]; - 5 + 5 * ((pxcor + pycor + 1) / 2 / nGrid) ] ; - This was for giving some gradient (to emphasize, visually, that this is a simulation) - Obsolete. Looks unnecessarily confusing.
    ask outpatches
    [ if any? neighbors with [member? self inpatches] [set pcolor black] ] ]
end


to update_pcolor
  let x1 map [[x] -> x_patch * x ] RGB-pip2
  let x2 map [[x] -> (1 - x_patch) * x ] RGB-pip1 ;  If 4 components vs. 3 components matter is involved, consider this code:  let x2 map [[x] -> (1 - x_initial) * x ] but-last RGB-pip1
  set pcolor (map + x1 x2)
end


to-report retrieve-simul_info_string
  report (word retrieve-num-date-time "_" Calculation-Type " " Enzyme-Pair-Type " " worldLength "um "
      "D_enz" D_enz " D_pip" D_pip " "
      "k_on-off-cat-Km " k_kon " " k_koff " " k_kcat " " k_Km " "
      "p_on-off-cat-Km-cat_unitA " k_kon " " k_koff " " k_kcat " " k_Km " " psol_kcat_per_A " "
      nGrid "Grid k_S-dependence-" k_S-dependence " p_Sdep-" p_S-dependence " initKPX" string-KPX )
end


to-report retrieve-num-date-time  ; current date in numerical format, yyyy-mm-dd
  let $dt substring date-and-time 16 27  ; "21-May-2013"
  let $dt2 substring date-and-time 0 8  ; 01:19
  report (word (substring $dt 7 11)           ; yyyy
    "-" (retrieve-month-num substring $dt 3 6)  ; mm
    "-" (substring $dt 0 2)           ; dd
    "_" (substring $dt2 0 2) "_" (substring $dt2 3 5) "_" (substring $dt2 6 8))
end


to-report retrieve-month-num [ #mon ]
  let $index 1 + position #mon
    ["Jan""Feb""Mar""Apr""May""Jun""Jul""Aug""Sep""Oct""Nov""Dec"]
  report substring (word (100 + $index)) 1 3  ; force 2-digit string
end
@#$#@#$#@
GRAPHICS-WINDOW
144
225
552
634
-1
-1
8.0
1
15
1
1
1
0
0
0
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
17
15
79
76
NIL
reset
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
19
215
82
276
NIL
go
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
866
336
1066
456
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
867
212
1066
332
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
"KIN" 1.0 2 -955883 true "" "if Calculation-Type = \"deterministic\" [plotxy time sum ([patch_k_density] of inpatches) * patchLength ^ 2]\nif Calculation-Type = \"stochastic\" [plotxy time count kinases]\n"
"PPT" 1.0 2 -8990512 true "" "if Calculation-Type = \"deterministic\" [plotxy time sum ([patch_p_density] of inpatches) * patchLength ^ 2]\nif Calculation-Type = \"stochastic\" [plotxy time count pptases]\n"

INPUTBOX
558
444
646
504
k_kon
0.5
1
0
Number

INPUTBOX
557
511
646
571
k_koff
1.0
1
0
Number

INPUTBOX
557
577
645
637
k_kcat
1.0
1
0
Number

INPUTBOX
653
442
742
502
p_kon
0.3
1
0
Number

INPUTBOX
655
508
744
568
p_koff
0.5
1
0
Number

INPUTBOX
653
576
745
636
p_kcat
3.0
1
0
Number

SWITCH
608
182
711
215
show_enz?
show_enz?
0
1
-1000

INPUTBOX
563
290
626
350
D_pip
2.0
1
0
Number

SLIDER
718
182
821
215
enz_size
enz_size
0
1
0.38
0.01
1
NIL
HORIZONTAL

SLIDER
195
114
480
147
world_pixel_length
world_pixel_length
100
700
400.0
50
1
NIL
HORIZONTAL

SLIDER
952
22
1044
55
log_timestep
log_timestep
-6.5
1
-2.0
0.1
1
NIL
HORIZONTAL

INPUTBOX
86
15
173
75
loss_threshold
0.0
1
0
Number

MONITOR
755
470
845
515
k/p (init.kin.adv)
k_kon * k_kcat / k_Km / (p_kon * p_kcat / p_Km)
6
1
11

PLOT
867
461
1066
581
Max Pon-patch
time
NIL
0.0
0.0
0.0
0.1
true
true
"" ""
PENS
"KIN" 1.0 2 -5298144 true "" "plotxy time max [k_Pon] of patches"
"PPT" 1.0 2 -15302303 true "" "plotxy time max [p_Pon] of patches"

CHOOSER
333
649
425
694
nGrid
nGrid
1 2 3 4 5 6 7 9 10 12 15 16 18 20 24 25 27 30 34 36 40 41 50 51 59 60 64 66 70 89 90 96 100 150 200 250 350 500
22

MONITOR
633
305
754
350
NIL
pip-diffuse-out-portion
7
1
11

SWITCH
492
113
582
146
wrap?
wrap?
1
1
-1000

SWITCH
148
650
320
683
save_timelapse_img?
save_timelapse_img?
0
1
-1000

SWITCH
23
651
141
684
record_vid?
record_vid?
0
1
-1000

SWITCH
867
176
1067
209
save_plots?
save_plots?
0
1
-1000

INPUTBOX
1163
22
1257
82
endtime
180.005
1
0
Number

INPUTBOX
340
15
605
107
input-geometry-fname
C:\\Users\\Neil\\Dropbox\\Research\\20200723\\confinements\\50-snail6.png
1
0
String

MONITOR
1062
23
1154
68
NIL
timestep
4
1
11

INPUTBOX
334
697
426
757
worldLength
30.0
1
0
Number

INPUTBOX
952
58
1046
118
timestep_input
0.0
1
0
Number

MONITOR
755
519
845
564
k/p (mac. adv)
k_kon * k_kcat / k_Km / k_koff / (p_kon * p_kcat / p_Km / p_koff)
6
1
11

SWITCH
1067
87
1264
120
Enzymes-touch-the-edge?
Enzymes-touch-the-edge?
0
1
-1000

PLOT
1080
347
1277
479
Kinase count - Edge vs non-Edge
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
"Edge" 1.0 0 -2674135 true "" "if any? edgeInPatches [plotxy time (count kinases-on edgeInPatches) / (count edgeInPatches)]"
"Non-Edge" 1.0 0 -13791810 true "" "if count non_Edge_Patches != 0 [plotxy time (count kinases-on non_Edge_Patches) / (count non_Edge_Patches)]"

CHOOSER
847
20
940
65
d_t-input-option
d_t-input-option
"exponent" "direct"
0

CHOOSER
181
16
330
61
geometry-setup
geometry-setup
"None" "Confinement"
1

SLIDER
22
689
141
722
vid_rec_intval
vid_rec_intval
1
300
200.0
10
1
NIL
HORIZONTAL

BUTTON
1067
123
1264
156
NIL
run_edge_test
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
148
687
320
720
tlapse_interval
tlapse_interval
0.5
150
5.0
0.5
1
s
HORIZONTAL

TEXTBOX
873
759
993
777
on-rate: (/x*um2*s)
12
0.0
1

SWITCH
868
710
1066
743
disallow-too-large-dx?
disallow-too-large-dx?
0
1
-1000

PLOT
868
584
1068
704
Max Min dx_patch
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
"max" 1.0 0 -16777216 true "" "plotxy time max [dx_patch] of inpatches "
"min" 1.0 0 -7500403 true "" "plotxy time min [dx_patch] of inpatches "

MONITOR
439
648
529
693
patchLength
worldLength / nGrid
7
1
11

SWITCH
185
178
314
211
display_time?
display_time?
1
1
-1000

INPUTBOX
553
709
643
769
k_Km
0.3
1
0
Number

INPUTBOX
652
709
746
769
p_Km
0.3
1
0
Number

INPUTBOX
652
643
745
703
psol_kcat_per_A
0.19
1
0
Number

INPUTBOX
564
223
625
283
D_enz
0.2
1
0
Number

CHOOSER
570
375
667
420
k_S-dependence
k_S-dependence
"linear" "MM"
1

BUTTON
22
727
141
760
NIL
save-recorded
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
122
150
172
210
lambda
2.0
1
0
Number

BUTTON
18
116
172
149
NIL
perturb-pips
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
19
151
116
211
perturb-amplitude
0.5
1
0
Number

CHOOSER
677
376
781
421
p_S-dependence
p_S-dependence
"linear" "MM"
1

BUTTON
182
69
331
102
NIL
change-input-geometry-file
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
631
237
754
282
NIL
enz-diffuse-out-portion
7
1
11

CHOOSER
615
74
807
119
Calculation-Type
Calculation-Type
"stochastic" "deterministic"
0

CHOOSER
616
125
787
170
Enzyme-Pair-Type
Enzyme-Pair-Type
"mem K - mem P" "mem K - sol P"
1

INPUTBOX
318
149
475
209
string-KPX
[0 0 0.5]
1
0
String

BUTTON
479
179
597
212
NIL
set-KPT-default
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
17
79
172
112
initial-fluctuation?
initial-fluctuation?
1
1
-1000

MONITOR
866
127
1038
172
NIL
mean [x_patch] of inpatches
6
1
11

PLOT
1080
680
1280
800
Max - Min dx_patch
NIL
NIL
0.0
0.0
0.0
0.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plotxy time max [dx_patch] of inpatches - min [dx_patch] of inpatches"

TEXTBOX
487
150
603
176
Initial patch density of \n[kinase, pptase, x]
10
0.0
1

PLOT
1080
503
1280
653
Density of Enzymes
NIL
NIL
0.0
0.0
0.0
0.0
true
true
"" ""
PENS
"dKIN" 1.0 0 -408670 true "" "if Calculation-Type = \"deterministic\" [plotxy time mean ([patch_k_density] of inpatches) / (count inpatches * patchLength ^ 2)]\nif Calculation-Type = \"stochastic\" [plotxy time count kinases / (count inpatches * patchLength ^ 2)]\n"
"dPPT" 1.0 0 -6759204 true "" "if Calculation-Type = \"deterministic\" [plotxy time mean ([patch_p_density] of inpatches) / (count inpatches * patchLength ^ 2)]\nif Calculation-Type = \"stochastic\" [plotxy time count pptases / (count inpatches * patchLength ^ 2)]\n"

SWITCH
334
764
527
797
ten-percent-height-test?
ten-percent-height-test?
1
1
-1000

MONITOR
21
286
78
331
q
2 * pi / lambda
5
1
11

MONITOR
22
340
79
385
NIL
time
3
1
11

SWITCH
146
724
318
757
snap-w-and-wo-enz?
snap-w-and-wo-enz?
0
1
-1000

TEXTBOX
1991
683
2113
893
Good working condition:\n\n50Snail6\nd_t 0.01\nworldLength 10\n\nmemP: \non 300-15\noff 10-5\ncat 2-2\nKm 10-0.1\n\nsolP:\nk_kon 420\nkcatP_per_uniA 2\n
11
0.0
1

TEXTBOX
1994
409
2144
619
Good working condition:\n\n50Snail6\nd_t 0.01\nworldLength 10\n\nmemP: \non 300-10\t\noff 10-5\ncat 1-1\nKm 10-0.1\n\nsolP:\nk_kon 420\nkcatP_per_uniA 1\n
11
0.0
1

TEXTBOX
1609
175
1759
385
Good working condition:\n\n50Snail6\nd_t 0.01\nworldLength 60\nmemP: \non 0.2-0.04\noff 1-1\ncat 1-20\nKm 0.2-0.2\n\nsolP:\nk_kon 0.1\nkcatP_per_uniA 0.04
11
0.0
1

TEXTBOX
1616
406
1756
616
Good working condition:\n\n50Snail6\nd_t 0.01\nworldLength 60\n\nmemP: \non 0.2-0.04\noff 1-1\ncat 1-10\nKm 0.3-0.3\n\nsolP:\nk_kon 0.2\nkcatP_per_uniA 0.09
11
0.0
1

PLOT
1081
171
1280
336
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

TEXTBOX
1437
183
1587
379
Good working condition:\n\n50Snail6\nd_t \nworldLength 10\nmemP: \non 300-160\noff 1-4\ncat 0.3-0.1\nKm 10-0.2\n\nsolP:\nk_kon 100\nkcatP_per_uniA 0.85
11
0.0
1

TEXTBOX
1781
782
1931
992
Good working [lINEAR] condition:\n\n50Snail6\nd_t \nworldLength 10\nmemP: \non 300-160\noff 1-3\ncat 0.3-0.1\nKm 10-0.2\n\nsolP:\nk_kon 100\nkcatP_per_uniA 0.8
11
0.0
1

TEXTBOX
1434
432
1584
628
Good working condition:\n\n50Snail6\nd_t \nworldLength 30\nmemP: \non 0.5-0.15\noff 1-1\ncat 1-10\nKm 0.3-0.3\n\nsolP:\nk_kon 0.5\nkcatP_per_uniA 0.2
11
0.0
1

SWITCH
0
429
102
462
setsavedir?
setsavedir?
1
1
-1000

SWITCH
0
468
137
501
simple-save-name?
simple-save-name?
1
1
-1000

TEXTBOX
1617
682
1767
892
$$$ Good working condition:\n\n50Snail6\nd_t 0.01\nworldLength 60\n\nmemP: \non 0.2-0.06\noff 1-1\ncat 1-5\nKm 0.3-0.3\n\nsolP:\nk_kon 0.1\nkcatP_per_uniA 0.04
11
0.0
1

TEXTBOX
1767
50
1917
246
** Similar to filling-up exp **\n\n50Snail6\nd_t \nworldLength 10\nmemP: \non 50-160\noff 1-4\ncat 1-0.5\nKm 10-0.2\n\nsolP:\nk_kon 50\nkcatP_per_uniA 1.29
11
0.0
1

TEXTBOX
1763
281
1913
491
**~~ Similar to filling-up exp **~~\n\n50Snail6\nd_t \nworldLength 10\nmemP: \non 100-85\noff 4-8\ncat 4-0.5\nKm 10-0.1\n\nsolP:\nk_kon 50\nkcatP_per_uniA 1.15
11
0.0
1

TEXTBOX
1774
539
1924
749
**~~ Similar to filling-up exp **~~\n\n50Snail6\nd_t \nworldLength 20\nmemP: \non 100-85\noff 4-8\ncat 4-0.5\nKm 10-0.1\n\nsolP:\nk_kon 50\nkcatP_per_uniA 1.15
11
0.0
1

TEXTBOX
1930
107
2080
317
@@@ Good working condition:\n\n50Snail6\nd_t 0.01\nworldLength 60\n\nmemP: \non 0.1-0.02\noff 1-0.5\ncat 1-5\nKm 0.3-0.3\n\nsolP:\nk_kon 0.1\nkcatP_per_uniA 0.04
11
0.0
1

TEXTBOX
1583
21
1733
91
10 1 1 10\n0.25\nThis is a pretty consistent solP condition\nwlen10
11
0.0
1

TEXTBOX
1305
649
1455
747
wlen30\n1 1 1-0.2\n0.35-0.2\n\n1 1 0.5-0.2\n.3 .5 1.5-0.2\n\n
11
0.0
1

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
    <enumeratedValueSet variable="Enzyme-Pair-Type">
      <value value="&quot;mem K - sol P&quot;"/>
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
