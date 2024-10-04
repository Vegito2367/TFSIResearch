<script>
	import "../app.css";
  // You can add logic here to handle form submissions, data, etc.
  function submitForm() {
    alert('Form submitted successfully!');
  }
</script>

<style lang="postcss">
  /* Basic styling for the page */
  :global(body) {
    font-family: Arial, sans-serif;
    margin: 0;
    padding: 0;
    background-color: #f4f4f4;
  }

  header, footer {
    background-color: #00509e;
    color: white;
    text-align: center;
    padding: 1rem;
  }

  .MainPage {
    padding: 2rem;
    text-align: center;
  }
  

  .hero {
    background-color: #e3f2fd;
    padding: 4rem 0;
    margin-bottom: 2rem;
  }

  .hero h1 {
    font-size: 2.5rem;
    color: #00509e;
		padding-bottom: 1.5rem;
		font-family: Verdana, Geneva, Tahoma, sans-serif;
  }

	h2{
		font-size: 1.25rem;
		color: #00509e;
		font-family: var(--font-mono);
	}

	.intro p {
		font-size: 1.25rem;
		color: #333;
		margin-bottom: 2rem;
		text-align: justify;
		font-family: var(--font-mono);
	}

  .intro, .highlights, .contact {
    margin-bottom: 2rem;
  }

  .button {
    background-color: #00509e;
    color: white;
    padding: 0.75rem 1.5rem;
    text-decoration: none;
    border-radius: 5px;
    font-weight: bold;
    transition: background-color 0.3s ease;
  }

  .button:hover {
    background-color: #00376d;
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
  }

</style>
<main class="MainPage mainPadding">
  <section class="hero">
    <h1>Analysis of Structural Properties of Anion</h1>
    <h2><i>Tej Gumaste, Fynnegan L. Cooper, James D. Blakemore</i></h2>
    <h2 style="padding-bottom: 3rem;"><b>University Of Kansas</b></h2>
    <a href="/Graphs" class="button">View Graphs</a>
    <a href="https://www.github.com/Vegito2367/TFSIResearch" target="_blank" class="button">Visit Github</a>
    <a href="https://blakemore.ku.edu/" class="button" target="_blank">Blakemore Lab</a>
  </section>

  <section class="intro">
    <h2>Project Overview</h2>
    <p>Chemical data science is an emerging field within the broader field of chemistry that has numerous applications of high relevance to a variety of academic and industrial pursuits. With quantum computing and artificial intelligence becoming more mainstream, it is of great interest for not only the data scientist, but the chemist as well, to take advantage of these technologies to streamline data processing, analyze large data sets, and reveal new chemical insights that would be otherwise hidden by the sheer amount and depth/complexity of available data. In this project, modern computing and fundamental chemical analysis have been combined in order to pursue new insights into the structural behavior of the weakly coordinating anion bis(trifluoromethylsulfonyl)imide, otherwise known as TFSI. Taking a data science approach, published solid-state crystal structures available in the Cambridge Structural Database (CSD)<sup>1,2</sup> including one or more TFSI species of interest were categorized and statistically analyzed using software built for this purpose. The goal of this project from a chemical perspective was to determine the structural characteristics displayed by TFSI as inferred from the structural data. The goal of this project from a data-science perspective was to develop a new software program using Python to parse Crystallographic Information Files (CIF) obtained from the CSD<sup>3</sup> into statistically relevant information that could be compared across the individual structural data sets. This research endeavor aims to highlight the applications of data science to an otherwise foreign area of research (structural chemistry) and outline the capabilities of data processing for future structural and chemical investigations.</p>
  </section>
  <section style=" text-align:left">
    <h2 style="color: black; padding-bottom:5px">References:</h2>
    <p style="color: black; text-align:left;padding-bottom:5px">1. C. R. Groom, I. J. Bruno, M. P. Lightfoot and S. C. Ward, Acta Cryst. B, 2016, 72, 171-179.</p>
    <p style="color: black;padding-bottom:5px">2. The analysis reported here was carried out using results from the 2024.2 CSD release, dated July 2024. This release of the CSD contains more than 1.3 million individual structures from which our TFSI-containing data were selected. More information about the CSD can be obtained from the Cambridge Crystallographic Data Centre. </p>
    <a href="https://www.ccdc.cam.ac.uk/structures/?gad_source=1&gclid=Cj0KCQjwpP63BhDYARIsAOQkATZJQgNzaPX-APCRoPGvpYzN-0-ZvIL6cKtLP2MV2ilXQvlmEx-9rYIaAnx8EALw_wcB" style="color: var(--color-theme-2);">3. Cambridge Crystallographic Data Center Website</a>
  </section>
</main>

<footer>
  <p>Made by Tej Gumaste - EECS and Chemistry @ KU</p>
</footer>


