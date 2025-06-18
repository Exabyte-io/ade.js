"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const JobContextMixin_1 = require("@mat3ra/code/dist/js/context/JobContextMixin");
const MaterialsContextMixin_1 = require("@mat3ra/code/dist/js/context/MaterialsContextMixin");
const MethodDataContextMixin_1 = require("@mat3ra/code/dist/js/context/MethodDataContextMixin");
const WorkflowContextMixin_1 = require("@mat3ra/code/dist/js/context/WorkflowContextMixin");
const made_1 = require("@mat3ra/made");
const MaterialContextMixin_1 = require("../../mixins/MaterialContextMixin");
const providers_1 = require("../../providers");
class VASPContextProvider extends providers_1.ExecutableContextProvider {
    constructor(config) {
        super(config);
        this._material = undefined;
        this._materials = [];
        this.initJobContextMixin();
        this.initMaterialsContextMixin();
        this.initMethodDataContextMixin();
        this.initWorkflowContextMixin();
        this.initMaterialContextMixin();
    }
    // eslint-disable-next-line class-methods-use-this
    buildVASPContext(material) {
        return {
            // TODO: figure out whether we need two separate POSCARS, maybe one is enough
            POSCAR: material.getAsPOSCAR(true, true),
            POSCAR_WITH_CONSTRAINTS: material.getAsPOSCAR(true),
        };
    }
    getDataPerMaterial() {
        if (!this.materials || this.materials.length <= 1) return {};
        return { perMaterial: this.materials.map((material) => this.buildVASPContext(material)) };
    }
    /*
     * @NOTE: Overriding getData makes this provider "stateless", ie. delivering data from scratch each time and not
     *        considering the content of `this.data`, and `this.isEdited` field(s).
     */
    getData() {
        // consider adjusting so that below values are read from PlanewaveDataManager
        // ECUTWFC;
        // ECUTRHO;
        return {
            ...this.buildVASPContext(this.material),
            ...this.getDataPerMaterial(),
        };
    }
}
VASPContextProvider.Material = made_1.Made.Material;
exports.default = VASPContextProvider;
(0, MaterialContextMixin_1.materialContextMixin)(VASPContextProvider.prototype);
(0, MaterialsContextMixin_1.materialsContextMixin)(VASPContextProvider.prototype);
(0, MethodDataContextMixin_1.methodDataContextMixin)(VASPContextProvider.prototype);
(0, WorkflowContextMixin_1.workflowContextMixin)(VASPContextProvider.prototype);
(0, JobContextMixin_1.jobContextMixin)(VASPContextProvider.prototype);
