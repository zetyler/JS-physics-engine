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

class Vec2 {
    constructor(x, y) {
        this.x = x;
        this.y = y;
    }

    addVec(v) {
        this.x += v.x;
        this.y += v.y;
    }

    addMult(v, s) {
        this.addVec(this.multiplied(v, s));
    }

    static multiplied(v, s) {
        return new Vec2(v.x * s, v.y * s);
    }

    static cross(a, b) {
        return a.x * b.y - a.y * b.x;
    }
}

class Body {
    constructor(shape, x, y) {
        this.position = new Vec2(x, y);
        this.velocity = new Vec2(0, 0);
        this.force = new Vec2(0, 0);
        this.orient = (Math.random() - 0.5) * Math.PI;
        this.angularVelocity = 0;
        this.torque = 0;
        this.staticFriction = 0.5;
        this.dynamicFriction = 0.3;
        this.restitution = 0.2;
        this.mass = 0;
        this.invMass = 0;
        this.inertia = 0;
        this.invInertia = 0;
        this.shape = shape;
        //shape.body = this;
        //shape.initialize();
    }

    applyForce(force) {
        this.force.addVec(force);
    }

    applyImpulse(impulse, contactVector) {
        this.velocity.addMult(impulse, this.invMass);
        this.angularVelocity += this.invInertia * Vec2.cross(contactVector, impulse);
    }

    setStatic() {
        this.mass = 0;
        this.invMass = 0;
        this.inertia = 0;
        this.invInertia = 0;
    }

    setOrient(radians) {
        orient = radians;
        //this.shape.setOrient(radians);
    }
}