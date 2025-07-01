function injectHTMLWithScripts(targetId, html) {
    const container = document.getElementById(targetId);
    container.innerHTML = html;

    const scripts = container.querySelectorAll('script');
    scripts.forEach(oldScript => {
        const newScript = document.createElement('script');
        newScript.type = oldScript.type || 'text/javascript';
        if (oldScript.src) {
            newScript.src = oldScript.src;
        } else {
            newScript.textContent = oldScript.textContent;
        }
        oldScript.replaceWith(newScript);
    });
}

document.getElementById('networkForm').addEventListener('submit', async function(e) {
    e.preventDefault();

    const genes = document.getElementById('genes').value;
    const tissues = Array.from(document.querySelectorAll('input[type="checkbox"]:checked')).map(cb => cb.value);
    const mode = document.getElementById('mode').value;
    const p_threshold = document.getElementById('p_threshold').value;

    if (!genes.trim()) {
        alert("Please enter at least one gene.");
        return;
    }

    document.getElementById('loading').style.display = 'block';
    document.getElementById('networkResult').innerHTML = '';

    const response = await fetch('/generate-network', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ genes, tissues, mode, p_threshold })
    })
    .then(response => response.json())
    .then(result => {
        document.getElementById('loading').style.display = 'none';
        //document.getElementById('networkResult').innerHTML = result.html_plot;
        injectHTMLWithScripts('networkResult', result.html_plot);
    })
    .catch(error => {
        console.error('Error:', error);
        alert("An error occurred.");
    });

});