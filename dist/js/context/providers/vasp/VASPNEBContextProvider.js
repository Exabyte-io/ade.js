"use strict";
var __importDefault =
    (this && this.__importDefault) ||
    function (mod) {
        return mod && mod.__esModule ? mod : { default: mod };
    };
Object.defineProperty(exports, "__esModule", { value: true });
const JobContextMixin_1 = require("@mat3ra/code/dist/js/context/JobContextMixin");
const MaterialsContextMixin_1 = require("@mat3ra/code/dist/js/context/MaterialsContextMixin");
const MethodDataContextMixin_1 = require("@mat3ra/code/dist/js/context/MethodDataContextMixin");
const WorkflowContextMixin_1 = require("@mat3ra/code/dist/js/context/WorkflowContextMixin");
const MaterialContextMixin_1 = require("../../mixins/MaterialContextMixin");
const MaterialsSetContextMixin_1 = require("../../mixins/MaterialsSetContextMixin");
const providers_1 = require("../../providers");
const VASPContextProvider_1 = __importDefault(require("./VASPContextProvider"));
class VASPNEBContextProvider extends providers_1.ExecutableContextProvider {
    constructor(config) {
        super(config);
        this._materials = [];
        this.initMaterialContextMixin();
        this.initMaterialsContextMixin();
        this.initMaterialsSetContextMixin();
        this.initMethodDataContextMixin();
        this.initWorkflowContextMixin();
        this.initJobContextMixin();
    }
    getData() {
        const sortedMaterials = this.sortMaterialsByIndexInSet(this.materials);
        const VASPContexts = sortedMaterials.map((material) => {
            const context = { ...this.config.context, material };
            const config = { ...this.config, context };
            return new VASPContextProvider_1.default(config).getData();
        });
        return {
            FIRST_IMAGE: VASPContexts[0].POSCAR_WITH_CONSTRAINTS,
            LAST_IMAGE: VASPContexts[VASPContexts.length - 1].POSCAR_WITH_CONSTRAINTS,
            INTERMEDIATE_IMAGES: VASPContexts.slice(1, VASPContexts.length - 1).map(
                (data) => data.POSCAR_WITH_CONSTRAINTS,
            ),
        };
    }
}
exports.default = VASPNEBContextProvider;
(0, MaterialContextMixin_1.materialContextMixin)(VASPNEBContextProvider.prototype);
(0, MaterialsContextMixin_1.materialsContextMixin)(VASPNEBContextProvider.prototype);
(0, MaterialsSetContextMixin_1.materialsSetContextMixin)(VASPNEBContextProvider.prototype);
(0, MethodDataContextMixin_1.methodDataContextMixin)(VASPNEBContextProvider.prototype);
(0, WorkflowContextMixin_1.workflowContextMixin)(VASPNEBContextProvider.prototype);
(0, JobContextMixin_1.jobContextMixin)(VASPNEBContextProvider.prototype);
