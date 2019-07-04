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
        this.set(x, y);
    }

    set(x, y) {
        this.x = x;
        this.y = y;
    }

    set(v) {
        this.x = v.x;
        this.y = v.y;
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

    static newArray(length) {
        let arr = [];
        for (let i = 0; i < length; ++i) {
            arr.push(new Vec2(0, 0));
        }
        return arr;
    }
}

class Shape {
    constructor(params) {
        this.type = params.type;
        this.radius = 0;
        //this.body;
        //this.u = new Mat2();
        if (this.type === 'c') {
            this.radius = params.radius;
        }
        else if (this.type === 'p') {
            this.vertices = Vec2.newArray(8);
            this.normals = Vec2.newArray(8);
            this.vertexCount = 0;
            if (params.hw) {
                vertexCount = 4;
                vertices[0].set(-params.hw, -params.hh);
                vertices[1].set(params.hw, -params.hh);
                vertices[2].set(params.hw, params.hh);
                vertices[3].set(-params.hw, params.hh);
                normals[0].set(0, -1);
                normals[1].set(1, 0);
                normals[2].set(0, 1);
                normals[3].set(-1, 0);
            }
            /*else {
                let verts = params.vertices;
                // Find the right most point on the hull
                let rightMost = 0;
                let highestX = verts[0].x;
                for (let i = 1; i < verts.length; ++i) {
                    let x = verts[i].x;
                    if (x > highestX) {
                        highestX = x;
                        rightMost = i;
                    }
                    // If same x, take lowest y
                    else if (x == highestX) {
                        if (verts[i].y < verts[rightMost].y) {
                            rightMost = i;
                        }
                    }
                }

                let hull = new Array(8);
                let outCount = 0;
                let indexHull = rightMost;

                for (;;) {
                    hull[outCount] = indexHull;

                    // Find next index that wraps around hull (find the next most counter-clockwise vertex)
                    let nextHullIndex = 0;
                    for (let i = 1; i < verts.length; ++i) {
                        // Skip if the same coord
                        if (nextHullIndex = indexHull) {
                            nextHullIndex = i;
                            continue;
                        }

                        // Cross every set of three unique verts
                    }
                }
            }*/
        }
    }
    
    clone() {
        if (this.type === 'c') {
            return new Shape({
                type: 'c',
                radius: this.radius
            });
        }
        else if (this.type === 'p') {
            let s = new Shape({
                type: 'p'
            });
            //s.u.set(this.u);
            for (let i = 0; i < this.vertexCount; ++i) {
                s.vertices[i].set(this.vertices[i]);
                s.normals[i].set(this.normals[i]);
            }
            s.vertexCount = this.vertexCount;
            return s;
        }
    }

    initialize() {}
    computeMass(density) {}
    setOrient(radians) {}
    getType() {}
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