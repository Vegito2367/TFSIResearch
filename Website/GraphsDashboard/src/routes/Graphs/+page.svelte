<script>
  // Add your script logic here
  import GraphCard from "../graphCard.svelte";
  let message = "η–N Metal Binding to TFSI in the solid state";
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

<main>
  <section class="hero">
    <p>{message}</p>
    <p class="subtitle">Results Dashboard</p>
  </section>
  <div>
    {#each metalToGraphMapping as metalUnit}
    <GraphCard
      title="S–N–S Angle vs. Bond Length avg for Structures containing {metalUnit.metal}"
      src={metalUnit.fileLocation}
    />
  {/each}
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
</style>
