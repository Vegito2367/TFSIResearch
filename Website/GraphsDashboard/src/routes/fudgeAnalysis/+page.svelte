<script>
  import hist from "$lib/histograms/SensitivityAnalysis.png"
</script>
<main style="margin: 15px;">
  <center><img src={hist} alt="Metal Presence"  style="margin-bottom:10px;height:450px;"></center>
  <p style="margin-left: 12%;margin-right: 12%;font-size: large;font-family:var(--font-mono); font-weight:600;">
    Histogram showing the relationship between the number of structures computed to feature a TFSI that is N-bound to a metal as a function of the f-value. In our chosen computational procedure, a given metal was considered to be bound to the N-atom of the TFSI core structure if the sum of the covalent radii of the N-atom of TFSI and the relevant metal (M) plus an empirical tolerance factor (denoted f) of 1 Å was greater than the raw distance between N and M. Considering this approach, the value of f determines the threshold at which M-N bonding appears to become significant. This histogram demonstrates how the f-value of 1.0 was selected, considering that it allows computational capture of nearly all the structures which appear to indeed feature TFSI that is N-bound to M. 
  </p>
</main>