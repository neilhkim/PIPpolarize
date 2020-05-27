extensions [
  vid
  csv
  palette
]


globals [
  ; Internally set Constants ;
  outerpcolor_ref
  pip1_color
  pip2_color
  kinase_color
  pptase_color
  seqKinase_color

  ; Calculated Constants ;
  timestep
  dist_kinase ; distance travelled by kinase per timestep
  dist_pptase ; distance travelled by pptase per timestep
  PIP_flux_frac ; how much PIP get dispersed to neighbors
  k_off_prob_sto
  p_off_prob_sto
  nInnerPatches
  nOuterPatches
  min_n_dominated_patch
  patchLength

  ; === ; Calcualted Constants --- Context : Confinement --- ;
  innerpatches
  outerpatches
  edgeInPatches
  non_Edge_Patches
  edgeOutPatches

  ; Trackers - ever changing ;
  avg_x
  time
  Sim_count
  cnt_k_win
  pKinaseWin
  seKinaseWin
  polarz_flag
  n_RxnsPolarized
  polarized_duration
  reactionTLength

  ; Internal options ;
  pptase_move?
  pptase_bind?
  pptase_unbind?
  bKnbP_flag

  ; Special calculated constant (not used often) ;
  area1
  area3
  area5
  area7
  area9

  avgx_C ; center
  stdx_C ; center
  avgx_NW ; northwest
  stdx_NW ; northwest
  avgx_SE ; southeast
  stdx_SE ; southeast

  sum_x_C ; x
  sum_x2_C ; x^2
  sum_x_NW ; x
  sum_x2_NW ; x^2
  sum_x_SE ; x
  sum_x2_SE ; x^2
]


breed [ pip1s pip1 ]
breed [ pip2s pip2 ]
breed [ kinases kinase ]
breed [ pptases pptase ]
breed [ clocks clock ]


turtles-own [
  orig_heading lr_collision? ud_collision? orig_dx orig_dy
  mydist
  seq?
  ;k_or_p - If you don't want kinase always overlapping on top of pptase, then try: assigning them to the same breed, but internally distinguishing between K and P.
]


patches-own [
  k_on_prob_sto
  p_on_prob_sto
  n_kin1
  n_ppt1
  x_patch
  new_xpatch ; this is a must-have
  real_neighbors
  n_neighbors
  dx_patch
]


to setup_world_from_input_file
  ifelse geometry_load_type = "None"
  ; If 1
  [ set innerpatches patches
    set outerpatches no-patches
    ask innerpatches [ set pcolor grey ]  ]
  ; Else 1
  [ ifelse geometry_load_type = "Confinement"
    ; If 2
    [ import-pcolors-rgb geometry_fname
      set innerpatches patches with [pcolor != [0 0 0]] ; loading from image files, syntax like "black" does not work. Have to use RGB.
      set outerpatches patches with [pcolor = [0 0 0]] ]  ; loading from image files, syntax like "black" does not work. Have to use RGB.
    ; Else 2
    [ ifelse geometry_load_type = "ROIs"
      ; If 3
      [ import-pcolors-rgb geometry_fname
        display   wait 0.5
        set area1 patches with [(item 0 pcolor) = 231 ]
        set area3 patches with [(item 0 pcolor) = 155 ]
        set area5 patches with [(item 0 pcolor) = 239 ]
        set area7 patches with [(item 0 pcolor) = 133 ]
        set area9 patches with [(item 0 pcolor) = 176 ]
        set innerpatches patches
        set outerpatches [] ]
      ; Else 3
      [ user-message "geometry_load_type selection logical error." ]
    ]
  ]
  set nInnerPatches count innerpatches
  ifelse any? outerPatches [ set nOuterPatches count outerpatches ]
  [ set nOuterPatches 0 ]
  ;--------2019-10-07-Edge enz count--------------;
  set edgeInPatches innerPatches with [n_neighbors < 8];
  set non_Edge_Patches innerPatches with [not member? self edgeInPatches]
  ;---------------------------------;
end


to reset
  clear-all
  clear-turtles
  clear-all-plots
  clear-output

  ; Functional set-ups ; --- ; Setup space ;
  resize-world 0 (nGrid - 1) 0 (nGrid - 1)
  set-patch-size world_pixel_length / nGrid
  ifelse wrap? [ __change-topology true  true  ]
                [ __change-topology false false ]
  setup_world_from_input_file
  if display_time?
  [  create-clocks 1 [ set label (precision time 1) setxy (min-pxcor + nGrid / 5) (max-pycor - nGrid / 20) set size 0 ]  ]

  ; Fix constants (internal) ;
  set outerpcolor_ref brown
  set pptase_move? true
  set pptase_bind? true
  set pptase_unbind? true
  set bKnbP_flag false
  set pip1_color [0 100 255] ;  set pip1_color extract-rgb blue
  set pip2_color [255 200 0]
  set kinase_color [255 150 0]
  set pptase_color [0 0 200]
  set seqKinase_color [255 50 0]
  ;  set kinase_color [255 0 0] - Prev: Red
  ;  set pptase_color [0 170 150] - Prev : Teal

  ; Calculate constants ;
  ; Calculate constants ; === ; Setup time
  RESET-TICKS
  ifelse timestepUses = "log_timestep"
  ; If 1
  [ set timestep 10 ^ log_timestep ]
  ; Else 1
  [ ifelse timestepUses = "timestep_input"
    ; If 2
    [ set timestep timestep_input ]
    ; Else 2
    [ user-message "timestep setup logical error."  ]  ]

  set k_off_prob_sto kinase_off_rate * timestep
  set p_off_prob_sto pptase_off_rate * timestep
  if k_off_prob_sto >= 1 [user-message "k_off_prob_sto > 1"]
  if p_off_prob_sto >= 1 [user-message "p_off_prob_sto > 1"]
  set patchLength worldLength / nGrid

  ; Calculate constants ;  === ; Setup diffusion properties
  set dist_kinase sqrt(4 * D_kinase * timestep) * (1 / worldLength * nGrid)
  set dist_pptase sqrt(4 * D_pptase * timestep) * (1 / worldLength * nGrid)
  set PIP_flux_frac 1 - exp(-4 * D_pip * timestep / (worldLength / nGrid) ^ 2)

  ; Set Trackers (ever changing) ;
  set Sim_count 0
  set cnt_k_win 0
  set avg_x x_initial

  set_neighbors_and_outColors
  ask innerpatches [ set x_patch x_initial ]
  ask innerpatches [ update_pcolor ]
end


to init_polarization_checker
  ; Initialize polarization detection parameters
  set polarz_flag false
  let bool_thisOnePolarized false
  set polarized_duration 0
  set min_n_dominated_patch nInnerPatches * min_pipfrac_for_dominance
end


to go
  setup_savedir_tlapseMsg_startvid
  set n_RxnsPolarized 0
  let cstr_go ""

  repeat N_repeat
  [ ask kinases [ die ]
    ask pptases [ die ]
    ; Initialize times
    set Sim_count Sim_count + 1
    set time 0
    let next_snapshot_time 0
    ; Initialize Polaziation Checker
    init_polarization_checker
    ; Initialized Patches
    ask innerPatches [ set x_patch x_initial ]
    ask innerPatches [ update_pcolor ]
    set avg_x x_initial

    ; Loop exit - when an experiment ends
    while [ (avg_x > loss_threshold) and (avg_x < 1 - loss_threshold) and (time <= endtime)  ]
    [
      ; Snapshot condition check
      if save_timelapse_img? and time >= next_snapshot_time [
        export-view (word commonsavestr "_tlapse n" Sim_count " t" (precision time 1) ".png")
        set next_snapshot_time next_snapshot_time + tlapse_interval
      ]
      if record_vid?
      [ if ticks mod vid_rec_intval = 0 [ vid:record-view ] ]

      if display_time?
      [ ask clocks [ set label (precision time 2) ] ]
      set time time + timestep;

      ; All the main functions are below:
      ifelse det_kinetics?
      ; If
      [ detConvert ]
      ; Else
      [ unbind
        enz_transform
        convert
        move
        bind ]
      ask innerPatches [ update_pcolor ]
      set avg_x mean [x_patch] of innerPatches
      check_polarization
      tick
    ]
    ; After an experiment, check the winner
    check_for_winner

    ; Clear non-accumulative plots
    set-current-plot "Max on-rates (patch)" clear-plot
    set-current-plot "mean_x of ROIs" clear-plot
    set-current-plot "Kinase count - Edge vs non-Edge" clear-plot

    ; Polarization result thus far (over multiple trials) - summarize
    if polarized_duration > min_pol_duration [ set n_RxnsPolarized n_RxnsPolarized + 1 ] ;    set pPolarized n_RxnsPolarized / Sim_count ;    set sePolarized sqrt(pPolarized * (1 - pPolarized) / Sim_count)

    ; Reaction duration
    set reactionTLength time

    ; Saving the snapshot of the last state.
    if Sim_count < 5 and save_timelapse_img?
    [ export-view (word commonsavestr "_tlapse n" Sim_count " t" (floor time) ".png") ]
    ; Set commonsavestr
    set cstr_go commonsavestr ; function calling (uses current time to generate saved file names)
  ]
  ; Finish-up Recording
  finishup_recording (cstr_go)
end


to bind
  ask innerpatches
  [ if pptase_bind?
    [ set p_on_prob_sto (pptase_on_rate * (1 - x_patch) * patchLength ^ 2 * timestep)
      if no_feedback? [ set p_on_prob_sto  (pptase_on_rate * patchLength ^ 2 * timestep) ]
      if p_on_prob_sto > 1 [    user-message (word "p_on_prob_sto: " p_on_prob_sto)  ]
      if random-float 1 < p_on_prob_sto
      [ sprout-pptases 1
        [ set mydist dist_pptase
          set lr_collision? false
          set ud_collision? false
          set seq? false
          ifelse show_enz?
          ; If
          [ set shape "circle"
            set size enz_size
            set color pptase_color
            setxy xcor + random-float 1 - 0.5 ycor + random-float 1 - 0.5 ]
          ; Else
          [ set size 0 ]  ]  ]  ]

    set k_on_prob_sto (kinase_on_rate * x_patch * patchLength ^ 2 * timestep)
    if no_feedback? [ set k_on_prob_sto (kinase_on_rate * patchLength ^ 2 * timestep) ]
    if k_on_prob_sto > 1 [   user-message (word "k_on_prob_sto: " k_on_prob_sto)    ]
    if random-float 1 < k_on_prob_sto
    [ sprout-kinases 1
      [ set mydist dist_kinase
        set lr_collision? false
        set ud_collision? false
        set seq? false
        ifelse show_enz?
        ; If
        [ set shape "circle"
          set size enz_size
          set color kinase_color
          setxy xcor + random-float 1 - 0.5 ycor + random-float 1 - 0.5 ]
        ; Else
        [ set size 0 ]  ]  ]  ]
end


to move
  ; Move PIPs ;
  if PIP_flux_frac > 1 [ user-message "diffusion fraction is more than 1" ]
  custom_diffuse (PIP_flux_frac) ; PIP_flux_frac is the function input

  ; Move Enzymes ;
  ifelse enz_touch_edge?
  ; If
  [ ; Move carelessly
    ask turtles
    [ rt random 360
      ifelse (patch-ahead mydist = nobody or not member? (patch-ahead mydist) innerPatches)
      ; If
      [ while [patch-ahead (mydist / 10) != nobody and member? (patch-ahead (mydist / 10)) innerPatches]
        [fd mydist / 10 ] ]
      ; Else
      [ fd mydist ]  ] ]
  ; Else
  [ ; Move carefully
    ask turtles
    [ right random 360
      set lr_collision? false
      set ud_collision? false

      if (member? (patch-ahead mydist) outerpatches) or (patch-ahead mydist = nobody)
      [ set orig_heading heading
        set orig_dx dx
        set orig_dy dy
        ; Test collision to left/right (for left, dx will be negative)
        set heading 90
        if (member? (patch-ahead (orig_dx * mydist)) outerpatches) or (patch-ahead (orig_dx * mydist) = nobody )
        [ set lr_collision? true ]
        ; Test collision to up/down (for down, dx will be negative)
        set heading 0
        if (member? (patch-ahead (orig_dy * mydist)) outerpatches) or (patch-ahead (orig_dy * mydist) = nobody)
        [ set ud_collision? true ]
        ; There is a special case where only when moving in the original heading, it will meet an obstacle. ; It's when running to a protruded corner. Then, depending on the incident angle, the reflection direction is determined.
        set heading orig_heading
        if lr_collision? = false and ud_collision? = false and (member? (patch-ahead mydist) outerpatches)
        [ let patchx [pxcor] of patch-ahead mydist
          let patchy [pycor] of patch-ahead mydist
          ifelse (patchy - 0.5 - ycor) / (patchx - 0.5 - xcor) > dy / dx
          ; If
          [ set ud_collision? true ]
          ; Else
          [ set lr_collision? true ]
        ]
        ; Apply the above determined reflection
        set heading orig_heading
        if lr_collision?
        [ set heading (- heading)
          set lr_collision? false ]
        if ud_collision?
        [ set heading (180 - heading)
          set ud_collision? false ]  ]

      forward mydist

      ; Now there is a special case where doing all of the above still makes the enzyme go over the boundary. This is when 2+ reflections is needed. In this case, I say the timestep is too big.
      if (member? patch-here outerpatches)
      [ set heading heading + 180
        pen-down
        forward mydist
        user-message "Reduce timescale. Model detects more than 1 reflection per timestep."
        pen-erase ]  ]  ]
end


to enz_transform
  if enz_transform_opt = "sequestration"
  [ ask kinases
    [ ifelse seq? = true
      ; If: it is already sequestered
      [ if random-float 1 < deseq_rate * timestep
        [ ; De-sequestration ;
          set mydist dist_kinase
          set seq? false
          if show_enz?
          [ set color kinase_color   set size enz_size set shape "circle"  ]   ]
      ]
      [ ; Else: it is not sequestered
        if random-float 1 < seq_rate * x_patch * timestep
        [ ; Sequesteration ;
          set mydist dist_kinase * seq_D_factor
          set seq? true
          if show_enz? [ set color seqKinase_color    set size enz_size * 1.0 set shape "circle" ]  ]  ]  ]  ]
end


to convert
  ask innerpatches
  [ let Jsum 0
    let PPT_contribution 0
    let KIN_contribution 0
    ; kinase: k1 (type 1) or k2 (type 2)
    set n_kin1 count kinases-here
    set n_ppt1 count pptases-here
    let n_kin2 count kinases-here with [seq? = true]

    let pptkcat 0
    ifelse bKnbP_flag
    ; If
;    [ set pptkcat ppt_kcat_patch ]
    [ set pptkcat ppt_kcat_per_um2 * patchLength ^ 2  ]
    ; Else
    [ set pptkcat pptase_kcat ]

    ifelse consider_Km?
    ; If
    [ set KIN_contribution   kinase_kcat * (1 - x_patch) / (kinase_Km + (1 - x_patch)) * n_kin1 / patchLength ^ 2                                             ; Contribution of the kinase type 1
      set KIN_contribution   KIN_contribution + (    kinase_kcat * kcatK_change * (1 - x_patch) / (kinase_Km + (1 - x_patch)) * n_kin2 / patchLength ^ 2    ) ; Contribution of the kinase type 2
      set PPT_contribution   -1 * pptkcat * x_patch / (pptase_Km + x_patch) * n_ppt1 / patchLength ^ 2  ]                                                 ; Contribution of the pptase type 1
    ; Else - Our current view (4/2/2020)
    [ set KIN_contribution   kinase_kcat * (1 - x_patch) * n_kin1 / patchLength ^ 2                                               ; Contribution of the kinase type 1
      set KIN_contribution   KIN_contribution + (    kinase_kcat * kcatK_change * (1 - x_patch) * n_kin2 / patchLength ^ 2      ) ; Contribution of the kinase type 2
      set PPT_contribution   -1 * pptkcat * x_patch * n_ppt1 / patchLength ^ 2   ]                                            ; Contribution of the pptase type 2

    set Jsum (KIN_contribution + PPT_contribution)
    set dx_patch Jsum * timestep
    if d_x_patch_restriction? and abs(dx_patch) > 0.5
    [ ; Inspect patch pxcor pycor
      user-message (word "abs(d_x_patch) > 0.5 patch: " dx_patch " at " pxcor " " pycor)  ]
    ifelse dx_patch > 0
    ; If 1
    [ ifelse dx_patch < (1 - x_patch)
      ; If 2-1
      [ set x_patch (x_patch + dx_patch) ]
      ; Else 2-1
      [ set x_patch 1  ] ]
    ; Else 1
    [ ; dx_patch <= 0
      ifelse -1 * dx_patch < x_patch
      ; If 2-2
      [  set x_patch (x_patch + dx_patch)  ]
      ; Else 2-2
      [  set x_patch 0 ]  ]  ]
end


to unbind
  ask kinases
  [ if seq? = false [
    if random-float 1 < k_off_prob_sto [ die ]  ]  ]
  if pptase_unbind?
  [ ask pptases
    [ if random-float 1 < p_off_prob_sto [ die ]  ]  ]
end

; 4/2/2020 이하로 보기.

to run_bound_kinase
  reset
  ; Fix unmoving phosphatases
  set pptase_move? false
  fix_1ppt_foreach_inpatch
  ; Randomly distribute N kinases
  put_n_kins_randomly (n_conj_kinase)
  setup_savedir_tlapseMsg_startvid
  while [time < endtime][
    if record_vid? and ticks mod vid_rec_intval = 0
    [ vid:record-view ]
    tick
    set time (time + timestep)
    convert
    move
    ask innerpatches [ update_pcolor ]
    set avg_x mean [x_patch] of innerpatches
    check_polarization
    update_ROI_avg_std
  ]
  if record_vid? [vid:save-recording (word commonsavestr "_conj_kin_vid.mp4") ]
  if save_interf_and_select_plots? [
    export-interface (word commonsavestr "_interface.png")
    export-plot "mean_x of ROIs" (word commonsavestr "_mean_x of ROIs.csv")
  ]
end


to setup_savedir_tlapseMsg_startvid
  if save_timelapse_img? or record_vid? or save_interf_and_select_plots?
  [ carefully [ set-current-directory user-directory ]
  ; Error-handling - do nothing at this time.
  [  ]  ]
  if save_timelapse_img?
  [ if N_repeat > 1
    [ user-message "Timelapse imaging for N_repeat > 1. Okay?" ]  ]
  if record_vid?
  [ vid:start-recorder ]
end

to finishup_recording [ inputstr ]
  if record_vid?     [
    vid:save-recording (word inputstr "_mov.mp4")
  ]
  if save_interf_and_select_plots? [
    export-interface (word inputstr "_interface.png")
    export-plot "Number of Enzymes" (word inputstr "_Number of Enzymes.csv")
    export-plot "PIP fraction" (word inputstr "_PIP fraction.csv")
  ]
end

to run_binding_K_vs_nonbinding_P
  reset
  set bKnbP_flag true
  fix_1ppt_foreach_inpatch
  set pptase_move? false
  set pptase_bind? false
  set pptase_unbind? false

  setup_savedir_tlapseMsg_startvid
  let cstr_bKIN_vs_nbPPT commonsavestr

  repeat N_repeat
  [ ask kinases [die]
    ; Initialize times
    set Sim_count Sim_count + 1
    set time 0
    let next_snapshot_time 0
    ; Initialized Patches
    ask innerPatches [ set x_patch x_initial ]
    ask innerPatches [ update_pcolor ]
    set avg_x x_initial

    while [time <= endtime]
    [
      if display_time?
      [ ask clocks [ set label (precision time 2) ] ]
      if save_timelapse_img? and time >= next_snapshot_time
      [ export-view (word cstr_bKIN_vs_nbPPT "_tlapse n" Sim_count " t" (precision time 1) ".png")
        set next_snapshot_time next_snapshot_time + tlapse_interval  ]
      ; Record video? - check
      if record_vid?
      [ if ticks mod vid_rec_intval = 0
        [ vid:record-view  ]  ]

      set time (time + timestep)
      ifelse det_kinetics?
      ; If
      [ det_convert_bK_vs_nbP ]
      ; Else
      [ unbind
        enz_transform
        convert
        move
        bind ]
      ask innerpatches [ update_pcolor ]
      set avg_x mean [x_patch] of innerpatches
      check_polarization
      tick
    ]
  ]
  finishup_recording (cstr_bKIN_vs_nbPPT)
end

to put_n_kins_randomly [ n_kin ]
  create-kinases n_kin
  [ set mydist dist_kinase
    set lr_collision? false
    set ud_collision? false
    set shape "circle"
    set size enz_size
    set color kinase_color
    setxy random-xcor random-ycor  ]
end

to fix_1ppt_foreach_inpatch
  ask innerpatches
  [ sprout-pptases 1
    [ set shape "circle"
      ifelse show_enz? and show_bPPT?
      ; If
      [ set size enz_size ]
      ; Else
      [ set size 0 ]
      set color pptase_color  ]  ]
end

to set_neighbors_and_outColors
  ask innerpatches
  [ set real_neighbors neighbors with [member? self innerPatches]
    set n_neighbors count real_neighbors  ]
  ; Decorate the outerpatches to make it visually obvious that this is simulation ;
  if geometry_load_type = "Confinement"
  [ ask outerpatches
    [ set pcolor outerpcolor_ref - 1 ]; - 5 + 5 * ((pxcor + pycor + 1) / 2 / nGrid) ] ; - This was for giving some gradient (to emphasize, visually, that this is a simulation) - Obsolete. Looks unnecessarily confusing.
    ask outerpatches
    [ if any? neighbors with [member? self innerPatches] [set pcolor black] ] ]
end


to custom_diffuse [ reducing_frac ]
  ask innerpatches
  [ set new_xpatch x_patch * (1 - reducing_frac * (n_neighbors / 8)) + (sum [ x_patch ] of real_neighbors) / 8 * reducing_frac  ]
  ask innerpatches
  [ set x_patch new_xpatch  ]
end

to check_polarization
  ifelse (count innerpatches with [(1 - x_patch) > min_pipfrac_for_dominance] / nInnerPatches > min_patchfrac_dominated ; threshold fraction of patches are dominated by PIP1
    and count innerpatches with [x_patch > min_pipfrac_for_dominance] / nInnerPatches > min_patchfrac_dominated)        ; threshold fraction of patches are dominated by PIP2
  ; If
  [ if polarz_flag = false
    [ set polarz_flag true
      set polarized_duration polarized_duration + timestep
      if polarized_duration > min_pol_duration
      [ output-print "Polarized" ] ] ]
  ; Else
  [ set polarz_flag false ]
end

to check_for_winner
    ; ==== IF pip2 is about to dominate === ;
    if avg_x > 1 - loss_threshold
    [ set cnt_k_win cnt_k_win + 1 ]
    set pKinaseWin cnt_k_win / Sim_count
    set seKinaseWin sqrt(pKinaseWin * (1 - pKinaseWin) / Sim_count)
end

to update_pcolor
  ifelse color_coding = "Blue and Yellow"
  ; If 1
  [ let x1 map [[x] -> x_patch * x ] pip2_color
    let x2 map [[x] -> (1 - x_patch) * x ] pip1_color ;  If 4 components vs. 3 components matter is involved, consider this:  let x2 map [[x] -> (1 - x_initial) * x ] but-last pip1_color
    set pcolor (map + x1 x2) ]
  ; Else 1
  [ ifelse color_coding = "Div - Red and Blue"
    ; If 2
    [ set pcolor palette:scale-gradient palette:scheme-colors "Divergent" "RdBu" 5 x_patch 1 0 ]
    ; Else 2
    [ user-message "color coding case error."]
  ]
end


to detConvert
  ask innerpatches
  [ let Jsum 0
    let PPT_contribution 0
    let KIN_contribution 0

    ifelse consider_Km?
    ; If
    [ set KIN_contribution   kinase_kcat * (1 - x_patch) / (kinase_Km + (1 - x_patch)) * (x_patch * kinase_on_rate / kinase_off_rate) / patchLength ^ 2
      set PPT_contribution   -1 * pptase_kcat * x_patch / (pptase_Km + x_patch) * ((1 - x_patch) * pptase_on_rate / pptase_off_rate) / patchLength ^ 2  ]
    ; Else
    [ set KIN_contribution   kinase_kcat * (1 - x_patch) * (x_patch * kinase_on_rate / kinase_off_rate) / patchLength ^ 2
      set PPT_contribution   -1 * pptase_kcat * x_patch * ((1 - x_patch) * pptase_on_rate / pptase_off_rate) / patchLength ^ 2  ]

    set Jsum    (KIN_contribution + PPT_contribution);
    set dx_patch Jsum * timestep
    if d_x_patch_restriction? and dx_patch > 0.5 [
      user-message (word "d_x_patch > 0.5 patch: " pxcor " " pycor)
    ]
    ifelse dx_patch > 0
    ; If 1
    [ ifelse dx_patch < (1 - x_patch)
      ; If 1-1
      [ set x_patch (x_patch + dx_patch)   ]
      ; Else 1-1
      [ set x_patch 1  ]  ]
    ; Else 1
    [ ; dx_patch <= 0
      ifelse -1 * dx_patch < x_patch
      ; If 1-2
      [ set x_patch (x_patch + dx_patch)  ]
      ; Else 1-2
      [ set x_patch 0  ]
    ]
  ]
end


to det_convert_bK_vs_nbP
  ask innerpatches
  [ let Jsum 0
    let PPT_contribution 0
    let KIN_contribution 0

    ifelse consider_Km?
    ; If
    [ set KIN_contribution   kinase_kcat * (1 - x_patch) / (kinase_Km + (1 - x_patch)) * (x_patch * kinase_on_rate / kinase_off_rate) / patchLength ^ 2
      set PPT_contribution   -1 * pptase_kcat * x_patch / (pptase_Km + x_patch) / patchLength ^ 2  ] ; <- PPTASE from solution
    ; Else
    [ set KIN_contribution   kinase_kcat * (1 - x_patch) * (x_patch * kinase_on_rate / kinase_off_rate) / patchLength ^ 2
      set PPT_contribution   -1 * pptase_kcat * x_patch / patchLength ^ 2  ]   ; ; <- PPTASE from solution

    set Jsum    (KIN_contribution + PPT_contribution);
    set dx_patch Jsum * timestep
    if d_x_patch_restriction? and dx_patch > 0.5
    [ user-message (word "d_x_patch > 0.5 patch: " pxcor " " pycor) ]

    ifelse dx_patch > 0
    ; If 1
    [ ifelse dx_patch < (1 - x_patch)
      ; If 1-1
      [ set x_patch (x_patch + dx_patch) ]
      ; Else 1-2
      [ set x_patch 1 ] ]
    ; Else 1
    [ ifelse -1 * dx_patch < x_patch
      ; If 1-2
      [ set x_patch (x_patch + dx_patch) ]
      ; Else 1-2
      [ set x_patch 0 ] ]
  ]
end

to run_edge_test
  reset
  if count edgeInPatches = 0
  [ user-message "There are no edge patches. Check 'wrap' setting." ]
  put_n_kins_randomly (n_conj_kinase)
  while [time < endtime]
  [ tick
    set time (time + timestep)
    move ]
end

to-report commonsavestr
  ifelse bKnbP_flag
  ; If 1
  [ ifelse consider_Km?
    ; If 1-2
    [ report (word num-date-time "_" worldLength "um "
      "D_ki " D_kinase ", D_pp " D_pptase " and D_pip " D_pip " "
      "k_on-off-cat-Km " kinase_on_rate " " kinase_off_rate " " kinase_kcat " " k_Km " "
      "p_cat_patch-Km " ppt_kcat_patch " " p_Km " "
      nGrid "Grid")
    ]
    ; Else 1-2
    [ report (word num-date-time "_" worldLength "um "
      "D_ki " D_kinase ", D_pp " D_pptase " and D_pip " D_pip " "
      "k_on-off-cat " kinase_on_rate " " kinase_off_rate " " kinase_kcat " "
      "p_cat_patch " ppt_kcat_patch " "
      nGrid "Grid") ] ]
  ; Else 1
  [
    ifelse consider_Km?
    ; If 1-1
    [ report (word num-date-time "_" worldLength "um "
      "D_ki " D_kinase ", D_pp " D_pptase " and D_pip " D_pip " "
      "k_on-off-cat-Km " kinase_on_rate " " kinase_off_rate " " kinase_kcat " " k_Km " "
      "p_on-off-cat-Km " pptase_on_rate " " pptase_off_rate " " pptase_kcat " " p_Km " "
      nGrid "Grid")
    ]
    ; Else 1-1
    [ report (word num-date-time "_" worldLength "um "
      "D_ki " D_kinase ", D_pp " D_pptase " and D_pip " D_pip " "
      "k_on-off-cat " kinase_on_rate " " kinase_off_rate " " kinase_kcat " "
      "p_on-off-cat " pptase_on_rate " " pptase_off_rate " " pptase_kcat " "
      nGrid "Grid") ] ]
end

to-report num-date-time  ; current date in numerical format, yyyy-mm-dd
  let $dt substring date-and-time 16 27  ; "21-May-2013"
  let $dt2 substring date-and-time 0 8  ; 01:19
  report (word (substring $dt 7 11)           ; yyyy
    "-" (month-num substring $dt 3 6)  ; mm
    "-" (substring $dt 0 2)           ; dd
    "_" (substring $dt2 0 2) "_" (substring $dt2 3 5) "_" (substring $dt2 6 8))
end

to-report month-num [ #mon ]
  let $index 1 + position #mon
    ["Jan""Feb""Mar""Apr""May""Jun""Jul""Aug""Sep""Oct""Nov""Dec"]
  report substring (word (100 + $index)) 1 3  ; force 2-digit string
end

to update_ROI_avg_std
  set sum_x_C      sum_x_C + mean [X_patch] of area5
  set sum_x2_C      sum_x2_C + (mean [X_patch] of area5) ^ 2
  set sum_x_NW      sum_x_NW + mean [X_patch] of area1
  set sum_x2_NW      sum_x2_NW + (mean [X_patch] of area1) ^ 2
  set sum_x_SE      sum_x_SE + mean [X_patch] of area9
  set sum_x2_SE      sum_x2_SE + (mean [X_patch] of area9) ^ 2

  set avgx_C     sum_x_C / ticks
  carefully [ set stdx_C     sqrt((sum_x2_C / ticks) - (avgx_C ^ 2)) ] [ set stdx_C 0 ]
  set avgx_NW     sum_x_NW / ticks
  carefully [ set stdx_NW     sqrt((sum_x2_NW / ticks) - (avgx_NW ^ 2)) ] [ set stdx_NW 0 ]
  set avgx_SE     sum_x_SE / ticks
  carefully [ set stdx_SE     sqrt((sum_x2_SE / ticks) - (avgx_SE ^ 2)) ] [ set stdx_SE 0 ]
end
@#$#@#$#@
GRAPHICS-WINDOW
75
197
483
606
-1
-1
20.0
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
19
0
19
1
1
1
ticks
30.0

BUTTON
55
10
117
43
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
54
44
117
77
NIL
go
NIL
1
T
OBSERVER
NIL
G
NIL
NIL
1

PLOT
548
390
748
510
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
"PIP1_1" 1.0 2 -13345367 true "" "if plots_opt = \"parallel\" [plotxy time 1 - avg_x]\nif plots_opt = \"sequential\" [plot 1 - avg_x]"
"PIP2_1" 1.0 2 -4079321 true "" "if plots_opt = \"parallel\" [plotxy time avg_x]\nif plots_opt = \"sequential\" [plot avg_x]"

MONITOR
752
121
816
166
NIL
time
4
1
11

PLOT
549
266
748
386
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
"KIN_1" 1.0 2 -955883 true "" "if plots_opt = \"parallel\" [plotxy time count kinases]\nif plots_opt = \"sequential\" [plot count kinases]"
"PPT_1" 1.0 2 -8990512 true "" "if plots_opt = \"parallel\" [plotxy time count pptases]\nif plots_opt = \"sequential\" [plot count pptases]"
"sKIN" 1.0 2 -2064490 true "" "if plots_opt = \"parallel\" [plotxy time count kinases with [seq? = true]]\nif plots_opt = \"sequential\" [plot count kinases with [seq? = true]]"

INPUTBOX
750
379
838
439
kinase_on_rate
3.0
1
0
Number

INPUTBOX
749
446
838
506
kinase_off_rate
2.0
1
0
Number

INPUTBOX
749
511
837
571
kinase_kcat
0.0
1
0
Number

INPUTBOX
845
377
934
437
pptase_on_rate
0.7
1
0
Number

INPUTBOX
847
443
936
503
pptase_off_rate
0.7
1
0
Number

INPUTBOX
845
511
937
571
pptase_kcat
4.5
1
0
Number

PLOT
1278
757
1438
877
pKinaseWin
Sim_count
NIL
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"pKinaseWin" 1.0 0 -16777216 true "" "plotxy Sim_count pKinaseWin"

MONITOR
1593
602
1664
647
NIL
pKinaseWin
4
1
11

SWITCH
923
74
1026
107
show_enz?
show_enz?
0
1
-1000

INPUTBOX
905
149
968
209
D_pip
2.0
1
0
Number

SLIDER
923
112
1026
145
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
56
107
341
140
world_pixel_length
world_pixel_length
100
700
400.0
100
1
NIL
HORIZONTAL

SLIDER
650
59
742
92
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
755
10
882
70
loss_threshold
0.001
1
0
Number

MONITOR
946
459
1036
504
k/p (init.kin.adv)
kinase_on_rate * kinase_kcat / kinase_Km / (pptase_on_rate * pptase_kcat / pptase_Km)
6
1
11

PLOT
549
515
748
635
Max on-rates (patch)
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
"KIN_1" 1.0 2 -5298144 true "" "if plots_opt = \"parallel\" [plotxy time max [k_on_prob_sto] of patches]\nif plots_opt = \"sequential\" [plot max [k_on_prob_sto] of patches]"
"PPT_1" 1.0 2 -15302303 true "" "if plots_opt = \"parallel\" [plotxy time max [p_on_prob_sto] of patches]\nif plots_opt = \"sequential\" [plot max [p_on_prob_sto] of patches]"

CHOOSER
546
12
638
57
nGrid
nGrid
1 2 3 4 5 6 7 9 10 12 15 16 18 20 24 25 27 30 34 36 40 41 50 51 59 66 89 90 96 100 150 200 250 350 500
13

CHOOSER
757
74
849
119
N_repeat
N_repeat
1 2 3 4 5 6 8 9 10 12 16 25 36 49 64 100 250 400 900 1000 1600
0

MONITOR
751
233
839
278
NIL
PIP_flux_frac
5
1
11

SWITCH
437
109
527
142
wrap?
wrap?
1
1
-1000

SWITCH
145
635
317
668
save_timelapse_img?
save_timelapse_img?
1
1
-1000

SWITCH
20
636
138
669
record_vid?
record_vid?
1
1
-1000

SWITCH
328
634
524
667
save_interf_and_select_plots?
save_interf_and_select_plots?
1
1
-1000

MONITOR
744
833
836
878
processivity_k
kinase_kcat / kinase_Km / kinase_off_rate
6
1
11

MONITOR
845
832
938
877
processivity_p
pptase_kcat / pptase_Km / pptase_off_rate
6
1
11

INPUTBOX
884
10
954
70
endtime
150.02
1
0
Number

INPUTBOX
264
11
529
103
geometry_fname
C:\\Users\\Neil\\OneDrive\\Netlogo - Shared Folder since 2019-10-22\\confinements\\250-snail6.png
1
0
String

MONITOR
649
164
741
209
NIL
timestep
4
1
11

MONITOR
1165
346
1257
391
Prob Polarized
n_RxnsPolarized / Sim_count
6
1
11

MONITOR
1166
395
1259
440
SE Polarized
sqrt((n_RxnsPolarized / Sim_count) * (1 - (n_RxnsPolarized / Sim_count)) / Sim_count)
5
1
11

MONITOR
1593
650
1664
695
seKinaseWin
seKinaseWin
4
1
11

INPUTBOX
547
60
639
120
worldLength
1.0
1
0
Number

MONITOR
1112
279
1216
324
NIL
polarized_duration
4
1
11

MONITOR
752
171
839
216
avg_x
avg_x
6
1
11

INPUTBOX
648
97
742
157
timestep_input
0.0
1
0
Number

PLOT
1466
728
1696
865
mean_x of ROIs
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
"🡔" 1.0 0 -1604481 true "" "if geometry_load_type = \"ROIs\" [plotxy time mean [x_patch] of area1]"
"🡕" 1.0 0 -6565750 true "" "if geometry_load_type = \"ROIs\" [plotxy time mean [x_patch] of area3]"
"C" 1.0 0 -16777216 true "" "if geometry_load_type = \"ROIs\" [plotxy time mean [x_patch] of area5]"
"🡗" 1.0 0 -11033397 true "" "if geometry_load_type = \"ROIs\" [plotxy time mean [x_patch] of area7]"
"🡖" 1.0 0 -6917194 true "" "if geometry_load_type = \"ROIs\" [plotxy time mean [x_patch] of area9]"

MONITOR
946
508
1036
553
k/p (mac. adv)
kinase_on_rate * kinase_kcat / kinase_Km / kinase_off_rate / (pptase_on_rate * pptase_kcat / pptase_Km / pptase_off_rate)
6
1
11

SWITCH
1125
23
1280
56
enz_touch_edge?
enz_touch_edge?
0
1
-1000

INPUTBOX
904
215
968
275
D_kinase
0.2
1
0
Number

INPUTBOX
971
214
1034
274
D_pptase
0.2
1
0
Number

PLOT
1114
74
1311
206
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
"Non-Edge_1" 1.0 0 -13791810 true "" "plotxy time (count kinases-on non_Edge_Patches) / (count non_Edge_Patches)"

CHOOSER
650
10
743
55
timestepUses
timestepUses
"log_timestep" "timestep_input"
0

CHOOSER
120
58
260
103
geometry_load_type
geometry_load_type
"None" "Confinement" "ROIs"
0

SLIDER
19
674
138
707
vid_rec_intval
vid_rec_intval
1
100
10.0
1
1
NIL
HORIZONTAL

CHOOSER
849
280
941
325
x_initial
x_initial
0 0.01 0.5 0.99 1
2

MONITOR
852
74
916
119
NIL
Sim_count
17
1
11

CHOOSER
1282
318
1374
363
n_conj_kinase
n_conj_kinase
1 10 100 1000 10000
0

BUTTON
1284
402
1399
436
NIL
run_bound_kinase
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1281
23
1360
56
NIL
run_edge_test
NIL
1
T
OBSERVER
NIL
E
NIL
NIL
1

BUTTON
1282
365
1451
399
NIL
run_binding_K_vs_nonbinding_P
NIL
1
T
OBSERVER
NIL
C
NIL
NIL
1

OUTPUT
1110
229
1215
271
11

SLIDER
1048
443
1232
476
min_pipfrac_for_dominance
min_pipfrac_for_dominance
0.5
1
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
1047
481
1231
514
min_patchfrac_dominated
min_patchfrac_dominated
0
1
0.1
0.1
1
NIL
HORIZONTAL

SLIDER
1047
519
1231
552
min_pol_duration
min_pol_duration
0
10
4.6
0.1
1
NIL
HORIZONTAL

SLIDER
145
672
317
705
tlapse_interval
tlapse_interval
0.5
150
1.0
0.5
1
s
HORIZONTAL

CHOOSER
121
10
261
55
color_coding
color_coding
"Blue and Yellow" "Div - Red and Blue"
0

CHOOSER
1411
12
1549
57
enz_transform_opt
enz_transform_opt
"default" "sequestration"
0

INPUTBOX
1418
68
1485
128
seq_rate
1.0
1
0
Number

TEXTBOX
1494
97
1525
115
/s/x
11
105.0
1

INPUTBOX
1418
135
1485
195
deseq_rate
0.5
1
0
Number

TEXTBOX
1492
160
1515
178
/s
11
105.0
0

INPUTBOX
1419
199
1484
259
seq_D_factor
0.5
1
0
Number

CHOOSER
549
219
687
264
plots_opt
plots_opt
"parallel" "sequential"
0

INPUTBOX
1420
266
1497
326
kcatK_change
1.0
1
0
Number

CHOOSER
775
333
913
378
rate_normalization
rate_normalization
"none" "normalize-by-x"
1

TEXTBOX
921
354
1041
372
on-rate: (/x*um2*s)
12
0.0
1

SWITCH
537
665
719
698
d_x_patch_restriction?
d_x_patch_restriction?
0
1
-1000

PLOT
1682
57
1882
207
Max Min dx_patch
NIL
NIL
0.0
10.0
0.0
0.0
true
false
"" ""
PENS
"max" 1.0 0 -16777216 true "" "plot max [dx_patch] of innerpatches"
"min" 1.0 0 -7500403 true "" "plot min [dx_patch] of innerpatches"

PLOT
1006
583
1206
733
Number of Enzymes (ind)
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"kin" 1.0 0 -955883 true "" "if plots_opt = \"parallel\" [plotxy time count kinases]\nif plots_opt = \"sequential\" [plot count kinases]"
"ppt" 1.0 0 -11221820 true "" "if plots_opt = \"parallel\" [plotxy time count pptases]\nif plots_opt = \"sequential\" [plot count pptases]"
"seqkin" 1.0 0 -2064490 true "" "if plots_opt = \"parallel\" [plotxy time count kinases with [seq? = true]]\nif plots_opt = \"sequential\" [plot count kinases with [seq? = true]]"

MONITOR
1273
466
1331
511
NIL
avgx_C
6
1
11

MONITOR
1338
467
1442
512
NIL
avgx_NW
6
1
11

MONITOR
1275
518
1333
563
NIL
stdx_C
6
1
11

MONITOR
1338
519
1441
564
NIL
stdx_NW
6
1
11

PLOT
1705
606
1906
726
Timed avgx
NIL
NIL
0.0
100.0
0.0
0.0
true
true
"" ""
PENS
"C" 1.0 0 -16777216 true "" "plot avgx_C"
"NW" 1.0 0 -2064490 true "" "plot avgx_NW"
"SE" 1.0 0 -8630108 true "" "plot avgx_SE"

PLOT
1706
728
1906
848
Timed stdev
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
"C" 1.0 0 -16777216 true "" "plot stdx_C"
"NW" 1.0 0 -2064490 true "" "plot stdx_NW"
"SE" 1.0 0 -8630108 true "" "plot stdx_SE"

TEXTBOX
379
686
529
728
5, 1.5, 3\n\n-,-,1.5
11
0.0
1

SWITCH
935
406
1048
439
no_feedback?
no_feedback?
1
1
-1000

TEXTBOX
785
686
894
714
if Km > 100,000 then V = kcat*n_Enz*[Sub]
11
0.0
1

MONITOR
548
124
638
169
patchLength
worldLength / nGrid
7
1
11

SWITCH
592
745
719
778
det_kinetics?
det_kinetics?
1
1
-1000

PLOT
1572
232
1772
382
reactionLength
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
"default" 1.0 0 -16777216 true "" "plotxy Sim_count reactionTLength"

INPUTBOX
426
751
498
811
kinase_Km
0.5
1
0
Number

INPUTBOX
505
752
572
812
pptase_Km
0.5
1
0
Number

SWITCH
777
788
913
821
consider_Km?
consider_Km?
1
1
-1000

INPUTBOX
1070
778
1164
838
ppt_kcat_patch
0.45
1
0
Number

SWITCH
60
144
189
177
display_time?
display_time?
0
1
-1000

INPUTBOX
764
718
830
778
k_Km
0.0
1
0
Number

INPUTBOX
865
719
925
779
p_Km
0.0
1
0
Number

INPUTBOX
840
577
945
637
ppt_kcat_per_um2
1.25
1
0
Number

SWITCH
973
37
1084
70
show_bPPT?
show_bPPT?
1
1
-1000

TEXTBOX
945
383
1095
401
was 0.7 (5/21/2020)
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
  <experiment name="experiment" repetitions="10" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>reset</setup>
    <go>go</go>
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
