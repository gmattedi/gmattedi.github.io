const savePdfButton = document.getElementById('btn-save-pdf');

if (savePdfButton) {
  savePdfButton.addEventListener('click', () => {
    window.print();
  });
}
