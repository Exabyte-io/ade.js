import _ from "underscore";
import lodash from "lodash";
import {mix} from "mixwith";
import s from "underscore.string";
import { Made } from "@exabyte-io/made.js";
import { PERIODIC_TABLE } from "@exabyte-io/periodic-table.js";

import {
    MaterialContextMixinBuilder,
    MethodDataContextMixin,
    WorkflowContextMixin,
    JobContextMixin
} from "@exabyte-io/code.js/dist/context";
import { ExecutableContextProvider } from "../../providers";

export class NWChemTotalEnergyContextProvider extends mix(ExecutableContextProvider).with(
    MaterialContextMixinBuilder(Made.Material),
    MethodDataContextMixin,
    WorkflowContextMixin,
    JobContextMixin,
) {

    get atomSymbols() {return this.material.Basis.uniqueElements}

    get atomicPositionsWithoutConstraints() {return this.material.Basis.atomicPositions}

    get atomicPositions() {return this.material.Basis.atomicPositionsWithConstraints};

    get cartesianAtomicPositions() {return this.material.toCartesian()};

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
        const BASIS = '6-31G';
        const FUNCTIONAL = 'B3LYP';

        return {
            CHARGE: CHARGE,
            MULT: MULT,
            BASIS: BASIS,
            NAT: this.atomicPositions.length,
            NTYP: this.atomSymbols.length,
            ATOMIC_POSITIONS: this.atomicPositions.join('\n'),
            ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS: this.atomicPositionsWithoutConstraints.join('\n'),
            ATOMIC_SPECIES: this.ATOMIC_SPECIES,
            FUNCTIONAL: FUNCTIONAL,
            CARTESIAN: this.cartesianAtomicPositions,
        }
    }

    get ATOMIC_SPECIES() {
        return _.map(this.atomSymbols, (symbol) => {
            return NWChemTotalEnergyContextProvider.symbolToAtomicSpecies(symbol);
        }).join('\n');
    }

    static symbolToAtomicSpecies(symbol, pseudo) {
        const el = PERIODIC_TABLE[symbol];
        const filename = pseudo ? lodash.get(pseudo, 'filename', s.strRightBack(pseudo.path, '/')) : '';
        return el ? s.sprintf('%s %f %s', symbol, el.atomic_mass, filename) : undefined;
    }
}