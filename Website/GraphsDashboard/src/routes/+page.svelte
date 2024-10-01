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

  main {
    padding: 2rem;
    text-align: center;
  }

  .hero {
    background-color: #e3f2fd;
    padding: 4rem 0;
    margin-bottom: 2rem;
  }

  .hero h1 {
    font-size: 3rem;
    color: #00509e;
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

  input, textarea {
    width: 100%;
    padding: 0.75rem;
    margin: 0.5rem 0;
    border: 1px solid #ccc;
    border-radius: 5px;
  }

  button {
    padding: 0.75rem 1.5rem;
    background-color: #00509e;
    color: white;
    border: none;
    border-radius: 5px;
    cursor: pointer;
  }

  button:hover {
    background-color: #00376d;
  }
</style>

<header>
  <h1>Chemistry Research Project</h1>
  <p>Exploring the molecular interactions in chemical reactions</p>
</header>

<main>
  <section class="hero">
    <h1>Welcome to Our Chemistry Research!</h1>
    <p>Join us as we delve into the molecular world and unravel the mysteries of chemical reactions.</p>
    <a href="#highlights" class="button">Explore Highlights</a>
  </section>

  <section class="intro">
    <h2>Project Overview</h2>
    <p>Our research focuses on understanding the molecular interactions that drive complex chemical reactions. With advanced tools and models, we aim to provide insights into reaction mechanisms that can impact various fields of chemistry, from industrial processes to environmental science.</p>
  </section>

  <section id="highlights" class="highlights">
    <h2>Research Highlights</h2>
    <div>
      <h3>Reaction Mechanism Analysis</h3>
      <p>We use computational tools to analyze the mechanisms behind key chemical reactions.</p>
    </div>
    <div>
      <h3>3D Modeling of Compounds</h3>
      <p>Our team employs state-of-the-art 3D modeling techniques to visualize molecular structures and interactions.</p>
    </div>
    <div>
      <h3>Collaboration with Industry</h3>
      <p>Partnering with industrial experts to apply our research in real-world chemical applications.</p>
    </div>
  </section>

  <section id="contact" class="contact">
    <h2>Contact Us</h2>
    <form on:submit|preventDefault={submitForm}>
      <input type="text" name="name" placeholder="Your Name" required />
      <input type="email" name="email" placeholder="Your Email" required />
      <textarea name="message" placeholder="Your Message" required></textarea>
      <button type="submit">Send Message</button>
    </form>
  </section>
</main>

<footer>
  <p>&copy; 2024 Chemistry Research Project | All Rights Reserved</p>
</footer>


