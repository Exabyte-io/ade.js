import { ContextProvider } from "@exabyte-io/code.js/dist/context";

export class ExecutableContextProvider extends ContextProvider {
    constructor(config) {
        super({
            ...config,
            domain: "executable"
        });
    }
}
