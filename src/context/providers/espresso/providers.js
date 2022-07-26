import _ from "underscore";
import lodash from "lodash";
import { mix } from "mixwith";
import s from "underscore.string";

import { Made } from "@exabyte-io/made.js";
import { PERIODIC_TABLE } from "@exabyte-io/periodic-table.js";

import {
    MaterialContextMixinBuilder,
    MaterialsContextMixinBuilder,
    MaterialsSetContextMixin,
    MethodDataContextMixin,
    WorkflowContextMixin,
    JobContextMixin
} from "@exabyte-io/code.js/dist/context";
import { ExecutableContextProvider } from "../../providers";

export class QEPWXContextProvider extends mix(ExecutableContextProvider).with(
    MaterialContextMixinBuilder(Made.Material),
    MethodDataContextMixin,
    WorkflowContextMixin,
    JobContextMixin,
) {

    get atomSymbols() {return this.material.Basis.uniqueElements}

    get atomicPositionsWithoutConstraints() {return this.material.Basis.atomicPositions}

    get atomicPositions() {return this.material.Basis.atomicPositionsWithConstraints};

    /*
     * @NOTE: Overriding getData makes this provider "stateless", ie. delivering data from scratch each time and not
     *        considering the content of `this.data`, and `this.isEdited` field(s).
     */
    getData() {

        // the below values are read from PlanewaveDataManager instead
        // ECUTWFC = 40;
        // ECUTRHO = 200;

        const IBRAV = 0;

        return {
            IBRAV,
            RESTART_MODE: this.RESTART_MODE,
            NAT: this.atomicPositions.length,
            NTYP: this.atomSymbols.length,
            ATOMIC_POSITIONS: this.atomicPositions.join('\n'),
            ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS: this.atomicPositionsWithoutConstraints.join('\n'),
            CELL_PARAMETERS: this.CELL_PARAMETERS,
            ATOMIC_SPECIES: this.ATOMIC_SPECIES,
        }
    }

    get RESTART_MODE() {
        return (this.job.parentJob || this.workflow.hasRelaxation) ? 'restart' : 'from_scratch';
    }

    getPseudoBySymbol(symbol) {
        return (this.methodData.pseudo || []).find(p => p.element === symbol);
    }

    get ATOMIC_SPECIES() {
        // atomic species with pseudopotentials
        return _.map(this.atomSymbols, (symbol) => {
            const pseudo = this.getPseudoBySymbol(symbol);
            return QEPWXContextProvider.symbolToAtomicSpecie(symbol, pseudo);
        }).join('\n');
    }

    get CELL_PARAMETERS() {
        return this.material.Lattice.vectorArrays.map(x => {
            return x.map(y => {
                return s.sprintf('%14.9f', y).trim();
            }).join(' ');
        }).join('\n');

    }

    static symbolToAtomicSpecie(symbol, pseudo) {
        const el = PERIODIC_TABLE[symbol];
        const filename = pseudo ? lodash.get(pseudo, 'filename', s.strRightBack(pseudo.path, '/')) : '';
        return el ? s.sprintf('%s %f %s', symbol, el.atomic_mass, filename) : undefined;
    }
}

export class QENEBContextProvider extends mix(ExecutableContextProvider).with(
    MaterialContextMixinBuilder(Made.Material),
    MaterialsContextMixinBuilder(Made.Material),
    MaterialsSetContextMixin,
    MethodDataContextMixin,
    WorkflowContextMixin,
    JobContextMixin,
) {

    getData() {
        const sortedMaterials = this.sortMaterialsByIndexInSet(this.materials);
        const PWXContexts = sortedMaterials.map(material => {
            const context = Object.assign({}, this.config.context, {material: material});
            const config = Object.assign({}, this.config, {context});
            return new QEPWXContextProvider(config).getData();
        });

        return {
            ..._.omit(PWXContexts[0], ["ATOMIC_POSITIONS", "ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS"]),
            FIRST_IMAGE: PWXContexts[0].ATOMIC_POSITIONS,
            LAST_IMAGE: PWXContexts[PWXContexts.length - 1].ATOMIC_POSITIONS,
            INTERMEDIATE_IMAGES: PWXContexts.slice(1, PWXContexts.length - 1).map(data => data.ATOMIC_POSITIONS),
        }
    }
}
