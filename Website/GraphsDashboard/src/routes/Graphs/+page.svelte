<script>
  // Add your script logic here
  import Graph from "../graph.svelte";
import GraphCard from "../graphCard.svelte";
  let message = "Analysis of Metal Binding to TFSI";
  let metalToGraphMapping = [
    {
      metal: "Au",
      fileLocation: "$lib/graphs/AuPres.html",
    },
    {
      metal: "Ag",
      fileLocation: "$lib/graphs/AgPres.html",
    },
    {
      metal: "Cu",
      fileLocation: "$lib/graphs/CuPres.html",
    },
    {
      metal: "Pt",
      fileLocation: "$lib/graphs/PtPres.html",
    },
    {
      metal: "Pd",
      fileLocation: "$lib/graphs/PdPres.html",
    },
    {
      metal: "Ru",
      fileLocation: "$lib/graphs/RuPres.html",
    },
    {
      metal: "Hg",
      fileLocation: "$lib/graphs/HgPres.html",
    },
    {
      metal: "Fe",
      fileLocation: "$lib/graphs/FePres.html",
    },
  ];
</script>

<main class="mainPadding">
  <section class="hero">
    <p>{message}</p>
    <p class="subtitle">Results Dashboard</p>
  </section>
  <div>
    {#each metalToGraphMapping as metalUnit}
    <GraphCard
      title="S–N–S Angle vs. Bond Length avg for Structures containing {metalUnit.metal}"
       element={metalUnit.metal}
    />
  {/each}
    <GraphCard
      title="Metal Presence in TFSI Structures"
      element="MetalPresenceHist" />

    <GraphCard
      title="TFSI structures metal-free vs metal containing" 
      element="MetalContainingHist" />
    
    <GraphCard
      title="Sensitivity Analysis"
      element="fudgeAnalysis" /> 
  </div>
  

</main>

<style>
  /* Add your styles here */
  p {
    color: #00509e;
    font-family: Verdana, Geneva, Tahoma, sans-serif;
    font-weight: 600;
    margin-bottom: 0;
    font-size: 2rem;
  }
  .hero {
    background-color: #e3f2fd;
  }

  div {
    display: flex;
    flex-wrap: wrap;
    justify-content: center;

  }
  
  .subtitle {
    font-size: 1.5rem;
    color: #00509e;
    font-family: Verdana, Geneva, Tahoma, sans-serif;
    font-weight: 300;
  }

  .mainPadding{
    flex: 1;
		display: flex;
		flex-direction: column;
		padding: 1rem;
		width: 100%;
		max-width: 64rem;
		margin: 0 auto;
		box-sizing: border-box;
    padding: 2rem;
    text-align: center;
  }
</style>
