import { ContextProvider } from "@exabyte-io/code.js/context";

export class ExecutableContextProvider extends ContextProvider {
    constructor(config) {
        super({
            ...config,
            domain: "executable",
        });
    }
}
