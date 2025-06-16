declare module "@exabyte-io/periodic-table.js" {
    interface PeriodicTableElement {
        atomic_mass: number;
        [key: string]: any;
    }

    export const PERIODIC_TABLE: Record<string, PeriodicTableElement>;
}
