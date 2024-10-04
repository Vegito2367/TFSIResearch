<script>
    // Add your script logic here
    export let element="Ag";
    let timestamp = Date.now();
</script>

<style>
    :global(body) {
        margin: 0;
        padding: 0;
    }

    .container {
        margin: 0;
        padding: 0;
    }
</style>

<div class="container">
    <!-- svelte-ignore a11y-missing-attribute -->
    <iframe src="https://vegito2367.github.io/HostedTFSIgraphs/{element}Graph.html?v={timestamp}" width="1500" height="900px"></iframe>
</div>