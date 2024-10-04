<script>
  // Add your script logic here
  export let title = "Graph Card";
  export let src = "https://example.com";
  export let element="Ag";
</script>

<style>
  /* Add your styles here */
  .card {
    border: 1px solid #ccc;
    border-radius: 8px;
    padding: 16px;
    box-shadow: 0 5px 7px rgba(0, 0, 0, 0.2);
    background-color: var(--offWhite);
    transition: transform 0.15s;
    margin: 10px;
    width: 15rem;
    height: 15rem;
  }

  .card:hover{
    box-shadow: 0 5px 7px rgba(0, 0, 0, 0.3);
    transform: translateY(-2px);
  }

  .title {
    font-size: 1.25em;
    margin-bottom: 15px;
    text-align: center;
    color: var(--color-theme-2);
    font-family: var(--font-mono);
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

</style>

<div class="card">
  <div class="title">{title}</div>
  <a href="/{element}" class="button" target="_blank">View Graph</a>
</div>