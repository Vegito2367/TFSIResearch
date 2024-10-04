<script>
 import Graph from "../graph.svelte";
 import { page } from "$app/stores";
 $:element = $page.params.element;
</script>
<main>
    <Graph element={element} />
</main>
