let canvas = document.querySelector('#render');
let ctx = canvas.getContext('2d');

function resize() {
    canvas.style.width = window.innerWidth + 'px';
    canvas.style.height = window.innerHeight + 'px';
    canvas.width = canvas.clientWidth;
    canvas.height = canvas.clientHeight;
}
window.addEventListener('resize', resize);
resize();

let keepRunning = true;
setTimeout(() => { 
    keepRunning = false; 
}, 10000);

let oldTime = (Date.now() / 1000);
function loop() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.lineWidth = 2;
    ctx.strokeStyle = 'hsl(200, 100%, 70%)';
    ctx.beginPath();
    ctx.moveTo(100, 100);
    ctx.lineTo(900, 600);
    ctx.closePath();
    ctx.stroke();

    ctx.font = "32px Arial";
    ctx.fillStyle = 'hsl(200, 100%, 70%)';
    let timeLeft = (10 - ((Date.now() / 1000) - oldTime)).toFixed(2);
    ctx.fillText(timeLeft, 50, 50);
    
    if(keepRunning)
        requestAnimationFrame(loop);
    else
        ctx.clearRect(0, 0, canvas.width, canvas.height);
}
loop();