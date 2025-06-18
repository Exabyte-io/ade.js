"use strict";
var __importDefault =
    (this && this.__importDefault) ||
    function (mod) {
        return mod && mod.__esModule ? mod : { default: mod };
    };
Object.defineProperty(exports, "__esModule", { value: true });
const periodic_table_js_1 = require("@exabyte-io/periodic-table.js");
const MethodDataContextMixin_1 = require("@mat3ra/code/dist/js/context/MethodDataContextMixin");
const made_1 = require("@mat3ra/made");
const lodash_1 = __importDefault(require("lodash"));
const underscore_1 = __importDefault(require("underscore"));
const underscore_string_1 = __importDefault(require("underscore.string"));
const JobContextMixin_1 = require("../../mixins/JobContextMixin");
const MaterialContextMixin_1 = require("../../mixins/MaterialContextMixin");
const WorkflowContextMixin_1 = require("../../mixins/WorkflowContextMixin");
const providers_1 = require("../../providers");
class NWChemTotalEnergyContextProvider extends providers_1.ExecutableContextProvider {
    constructor(config) {
        super(config);
        this._material = undefined;
        this.initMethodDataContextMixin();
        this.initWorkflowContextMixin();
        this.initJobContextMixin();
        this.initMaterialContextMixin();
    }
    get atomicPositionsWithoutConstraints() {
        return this.material.Basis.atomicPositions;
    }
    get atomicPositions() {
        return this.material.Basis.getAtomicPositionsWithConstraintsAsStrings();
    }
    get atomSymbols() {
        return this.material.Basis.uniqueElements;
    }
    get cartesianAtomicPositions() {
        return this.material.Basis.toCartesian !== undefined;
    }
    get ATOMIC_SPECIES() {
        return underscore_1.default
            .map(this.atomSymbols, (symbol) => {
                return NWChemTotalEnergyContextProvider.symbolToAtomicSpecies(symbol);
            })
            .join("\n");
    }
    /*
     * @NOTE: Overriding getData makes this provider "stateless", ie. delivering data from scratch each time and not
     *        considering the content of `this.data`, and `this.isEdited` field(s).
     */
    getData() {
        /*
        TODO: Create ability for user to define CHARGE, MULT, BASIS and FUNCTIONAL parameters.
         */
        const CHARGE = 0;
        const MULT = 1;
        const BASIS = "6-31G";
        const FUNCTIONAL = "B3LYP";
        return {
            CHARGE,
            MULT,
            BASIS,
            NAT: this.atomicPositions.length,
            NTYP: this.atomSymbols.length,
            ATOMIC_POSITIONS: this.atomicPositions.join("\n"),
            ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS: this.atomicPositionsWithoutConstraints.join("\n"),
            ATOMIC_SPECIES: this.ATOMIC_SPECIES,
            FUNCTIONAL,
            CARTESIAN: this.cartesianAtomicPositions,
        };
    }
    static symbolToAtomicSpecies(symbol, pseudo) {
        const el = periodic_table_js_1.PERIODIC_TABLE[symbol];
        const filename = pseudo
            ? lodash_1.default.get(
                  pseudo,
                  "filename",
                  underscore_string_1.default.strRightBack(pseudo.path || "", "/"),
              )
            : "";
        return el
            ? underscore_string_1.default.sprintf("%s %f %s", symbol, el.atomic_mass, filename)
            : undefined;
    }
}
NWChemTotalEnergyContextProvider.Material = made_1.Made.Material;
exports.default = NWChemTotalEnergyContextProvider;
(0, MaterialContextMixin_1.materialContextMixin)(NWChemTotalEnergyContextProvider.prototype);
(0, MethodDataContextMixin_1.methodDataContextMixin)(NWChemTotalEnergyContextProvider.prototype);
(0, WorkflowContextMixin_1.workflowContextMixin)(NWChemTotalEnergyContextProvider.prototype);
(0, JobContextMixin_1.jobContextMixin)(NWChemTotalEnergyContextProvider.prototype);
